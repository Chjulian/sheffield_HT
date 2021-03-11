## Identify a starting tree with the fewest number of 'imports' as given by the
## timed contact data

get_initial_tree <- function(data, config,
                             priors = custom_priors(),
                             likelihoods = custom_likelihoods(),
                             n_iter = 5e3,
                             max_dist = 5,
                             days_after_admission = 0,
                             scale = 1) {

  ## set up config to only look for direct ward contacts
  tmp_config <- create_config(
    n_iter = n_iter, init_pi = 1, move_pi = FALSE,
    init_eps = ifelse(config$init_eps == 1, 1, 0.9), move_eps = FALSE,
    n_iter_import = 2e3,
    init_tree = 'star', sample_every = 100,
    find_import = FALSE,
    move_model = FALSE, swap_place = TRUE,
    move_alpha = TRUE, move_t_inf = TRUE,
    move_joint = FALSE,
    init_kappa = rep(1, length(data$dates)),
    move_kappa = FALSE, move_tau = FALSE,
    p_wrong = 1e-15
  )
  data$p_wrong <- tmp_config$p_wrong

  ## do first run
  cat("Finding direct ward contacts...\n")
  first_run <- outbreaker(data, tmp_config, priors = priors, likelihoods = likelihoods)

  ## get false negatives
  fn <- get_ward_fn(first_run)

  ## carry over ancestries and infection times to second run
  n <- nrow(first_run)
  nm <- names(first_run)
  tmp_config$init_alpha <- tmp_config$init_tree <- unlist(first_run[n, grep('alpha', nm)])
  tmp_config$init_t_inf <- unlist(first_run[n, grep('t_inf', nm)])

  ## if a false negative, set kappa to 2
  tmp_config$init_kappa <- as.integer(ifelse(fn, 2, 1))
  tmp_config$init_t_onw <- tmp_config$init_t_inf - 1

  ## find those with only 1 day between t_inf and t_inf_alpha - this makes an
  ## intermediate case impossible so we set the starting infector to earliest
  ## case
  to_fix <- tmp_config$init_t_onw - tmp_config$init_t_inf[tmp_config$init_tree] == 0
  tmp_config$init_tree[which(to_fix)] <- which.min(data$dates)
  tmp_config$init_alpha <- tmp_config$init_tree

  ## keep the number of unobserved cases fixed but explore ancestries
  tmp_config[c("move_kappa", "move_alpha", "find_import", "move_tau", "swap_place")] <- FALSE
  tmp_config$move_joint <- tmp_config$move_model <- TRUE
  data$swap_place <- FALSE
  tmp_config$sd_t_onw <- 30*scale

  ## edge cases if tau or eps is 1
  tmp_config$init_pi <- 0.5
  if(config$init_tau == 1) {
    tmp_config$init_tau <- 1 - 1e-15
    data$p_wrong <- 1e-15
  } else {
    tmp_config$init_tau <- 0.5
  }
  if(config$init_eps == 1) {
    tmp_config$init_eps <- 1 - 1e-15
    data$p_wrong <- 1e-15
  } else {
    tmp_config$init_eps <- 0.5
  }

  ## second run
  cat("Finding indirect ward contacts...\n")
  sec_run <- outbreaker(data, tmp_config, priors = priors, likelihoods = likelihoods)
  if(config$init_eps == 1) sec_run$eps <- 1
  if(config$init_tau == 1) sec_run$tau <- 1

  ## get remaining false negatives
  fn <- get_ward_fn(sec_run) == 1

  ## cases with onset < x days after admission are also labelled as imports
  after_admission <- data$ctd_timed %>%
    group_by(id) %>%
    summarise(adm = min(adm)) %>%
    left_join(tibble(id = data$ids, onset = data$dates), ., "id") %>%
    ## if ward data is missing, set as import
    mutate(import = replace_na(onset - adm <= days_after_admission*scale, TRUE))

  ## set as imports
  fn[after_admission$import] <- TRUE

  ## carry over ancestries, infection times and kappa
  n <- nrow(sec_run)
  nm <- names(sec_run)
  config$init_alpha <- config$init_tree <- unlist(sec_run[n, grep('alpha', nm)])
  config$init_t_inf <- unlist(sec_run[n, grep('t_inf', nm)])
  config$init_t_onw <- unlist(sec_run[n, grep('t_onw',nm)])
  config$init_kappa <- unlist(sec_run[n, grep('kappa',nm)])
  config$init_t_onw[config$init_kappa == 1] <- -1000

  ## include max_dist in false negatives/imports
  fn[data$D[cbind(seq_along(config$init_alpha), config$init_alpha)] > max_dist] <- TRUE

  ## no longer find imports as they have been identified above
  config$find_import <- FALSE

  ## settings for false negatives
  config$init_alpha[fn] <-
    config$init_tree[fn] <-
    config$init_t_onw[fn] <-
    config$init_kappa[fn] <- NA

  cat(paste0("Identified ", sum(is.na(config$init_alpha)),
             " seperate transmission chains.\n"))

  ## first_run <<- first_run
  ## sec_run <<- sec_run
  ## fn <<- fn

  return(config)

}

## extract ward info
get_ward <- function(i, t, data) {
  ind <- t + data$C_ind + 1
  ind <- ind[ind %in% seq_len(ncol(data$ctd_timed_matrix[[1]]))]
  data$ctd_timed_matrix[[1]][i, ind]
}

## Get the number of false negatives (ie missing contacts) in your proposed tTree
get_ward_fn <- function(x) {

  calc_fn <- function(from, to, t_inf, t_onw, kappa, eps, tau) {

    if(is.na(from)) return(0)

    if(kappa == 1) {
      ind1 <- t_inf + C_ind + 1
      if(ind1 < 1 | ind1 > ncol(data$ctd_timed_matrix)) return(1)
      ## If eps == 1 then direct trasmission events have to be on the same ward
      if(eps == 1) {
        if(data$ctd_timed_matrix[from, ind1] == data$ctd_timed_matrix[to, ind1] &
           data$ctd_timed_matrix[to, ind1] != -1) return(0) else return(1)
      } else {
        if(data$ctd_timed_matrix[from, ind1] != -1 &
           data$ctd_timed_matrix[to, ind1] != -1) return(0) else return(1)
      }
    } else if(kappa > 1) {
      ind1 <- t_onw + C_ind + 1
      ind2 <- t_inf + C_ind + 1
      if(ind1 < 1 | ind1 > ncol(data$ctd_timed_matrix)) return(1)
      if(ind2 < 1 | ind1 > ncol(data$ctd_timed_matrix)) return(1)
      ## If tau == eps == 1, the infector, intermediate and infectee have to be
      ## on the same ward
      if(tau == 1 & eps == 1) {
        if(data$ctd_timed_matrix[from, ind1] == data$ctd_timed_matrix[to, ind2] &
           data$ctd_timed_matrix[to, ind2] != -1) return(0) else return(1)
      } else {
        if(data$ctd_timed_matrix[from, ind1] != -1 &
           data$ctd_timed_matrix[to, ind2] != -1) return(0) else return(1)
      }
    }
  }

  ## get dimensions and name
  n <- nrow(x)
  nm <- names(x)

  ## extract data
  from <- unlist(x[n, grep('alpha', nm)])
  to <- seq_along(from)
  t_inf <- unlist(x[n, grep('t_inf', nm)])
  t_onw <- unlist(x[n, grep('t_onw', nm)])
  kappa <- unlist(x[n, grep('kappa', nm)])
  eps <- unlist(x[n, grep('eps', nm)])
  tau <- unlist(x[n, grep('tau', nm)])

  if(is.null(eps)) eps <- 0.5
  if(is.null(tau)) tau <- 0.5
  if(is.null(t_onw)) t_onw <- rep(1, length(from))
  t_onw[is.na(t_onw)] <- -1

  data$ctd_timed_matrix <- data$ctd_timed_matrix[[1]]
  C_ind <- data$C_ind

  ## calculate false negatives
  out <- mapply(calc_fn, from, to, t_inf, t_onw, kappa, MoreArgs = list(eps, tau))

  return(out)

}

## Visualise the wards
vis_ward <- function(res, meta, wards, ord = NULL, consensus = FALSE,
                     alph = 0.01, size = 4, nrow = 2, thin = 1) {

  ind <- seq(1, nrow(res), by = thin)
  res <- res[ind,]

  labs <- as.character(meta$id)

  if(is.null(ord)) {
    ord <- get_order(res)
    tree_sort <- match(ord, labs)
    ord <- match(meta$id, ord)
  } else {
    tree_sort <- match(ord, labs)
    ord <- match(meta$id, ord)
  }

  lett <- apply(combn(letters[1:26], 2), 2, paste, collapse = "")
  meta$id <- lett[ord]
  wards$id <- meta$id[match(wards$id, labs)]

  ## Get inferred infection times
  t_inf <- as.matrix(res[,grep("t_inf", names(res))])
  alpha <- as.matrix(res[,grep("alpha", names(res))])
  t_onw <- as.matrix(res[,grep("t_onw", names(res))])
  kappa <- as.matrix(res[,grep("kappa", names(res))])

  id <- meta$id[col(t_inf)]
  t_inf <- as.vector(t_inf)
  t_onw <- as.vector(t_onw)
  kappa <- as.vector(kappa)
  alpha  <- meta$id[as.vector(alpha)]

  meta$id %<>% factor(levels = .[ord])
  wards$id %<>% factor(levels = levels(meta$id))
  id %<>% factor(levels = levels(meta$id))
  alpha %<>% factor(levels = levels(meta$id))

  if(!consensus) {
    df2 <- data.frame(date_from = t_inf,
                      date_to = t_inf,
                      kappa = kappa,
                      id = id,
                      alpha = alpha,
                      support = 1)
    my_scale <- scale_size_continuous(range = c(0.5, 0.5))
  } else {
    tree <- summary(res)$tree
    id <- meta$id
    alpha <- meta$id[tree$from]
    id %<>% factor(levels = levels(meta$id))
    alpha %<>% factor(levels = levels(meta$id))
    t_onw <- tree$onw
    df2 <- data.frame(date_from = tree$time,
                      date_to = tree$time,
                      kappa = tree$generations,
                      id = id,
                      alpha = alpha,
                      support = tree$support)
    my_scale <- scale_size_continuous(range = c(0.2, 1))
    alph <- 1
  }

  df2$date_from[df2$kappa > 1 & !is.na(df2$kappa)] <- t_onw[df2$kappa > 1 & !is.na(df2$kappa)]
  df2 %<>%
    mutate(
      date_from = as.Date(.$date_from, origin = min(meta$date_onset)),
      date_to = as.Date(.$date_to, origin = min(meta$date_onset)),
      kappa = .$kappa > 1
    )

  df2 <- df2 %>%
    group_by(date_from, date_to, kappa, id, alpha) %>%
    mutate(support = n()) %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(support = support/sum(support))


  df3 <- transmute(meta, id = id, date = date_onset)

  ## This will ensure ward occupancy is displayed until the end of day of discharge
  wards$dis %<>% add(1)

  ggplot(wards) +
    geom_segment(aes(y = adm - 0.5, x = id,
                     xend = id, yend = dis + 0.5),
                 color = 'black',
                 size = size + 0.5) +
    geom_segment(aes(y = adm, x = id, xend = id, yend = dis, colour = ward),
                 size = size) +
    geom_segment(data = df2, aes(x = alpha, xend = id, size = support,
                                 y = date_from, yend = date_to, linetype = kappa, alpha = support),
                 arrow = arrow(length = unit(0.1,"cm"))) +
    geom_point(data = df3, aes(id, date)) +
    coord_flip() +
    scale_x_discrete(labels = labs[tree_sort]) +
    scale_alpha_continuous(range = c(0.05, 0.5)) +
    scale_y_date(date_labels = "%b %Y",
                 limits = c(min(wards$adm - 0.5), NA)) +
    my_scale +
    labs(x = 'Case ID', y = 'Date') +
    theme_minimal(base_size = 11) +
    guides(linetype = FALSE,
           size = FALSE,
           color = guide_legend(nrow = nrow)) +
    scale_linetype_manual(values = c("solid", "solid")) +
    theme(legend.position = 'bottom', legend.direction = 'horizontal')

}

## Recursive function to identify how many layers deep a case is
get_depth <- function(i, counter, alpha) {
  ances <- alpha[i]
  if(is.na(ances)) {
    return(counter)
  } else {
    counter <- counter + 1
    get_depth(ances, counter, alpha)
  }
}

## This function will order the cases in such a manner that the infector is
## closest to its infectees, starting at the leaves of the tree and moving to
## the root
get_order <- function(res) {

  get_sq <- function(len) {
    if(len %% 2 != 0) {
      len <- len + 1
      odd <- TRUE
    } else {
      odd <- FALSE
    }
    sq <- seq(1, 1 + 2*(len - 1), 2)
    out <- sq - median(sq)
    if(odd) out <- out[-1]
    return(out)
  }

  tree <- summary(res)$tree
  epi <- make_epicontacts(tree, tree, id = 2, na_rm_contacts = TRUE)

  coor <- get_coor(
    epi, 't_inf',
    rank_contact = 'from',
    axis_type = 'none',
    unlinked_pos = 'bottom'
  )

  return(meta$id[order(coor$y)])

}

## summarise ward data (all from-to and ward_from-ward_to combinations)
get_ward_summary <- function(res, data) {

  ## extract information for a single MCMC draw
  extract <- function(from, to, t_inf, t_onw, kappa) {

    if(is.na(from)) {
      out <- tibble(from = from, to = to, ward_from = NA, ward_to = NA, kappa = NA, t_inf = NA)
    } else {
      if(kappa == 1) {
        ind1 <- t_inf + C_ind + 1
        ward_from <- data$ctd_timed_matrix[from, ind1]
        ward_to <- data$ctd_timed_matrix[to, ind1]
      } else if(kappa > 1) {
        ind1 <- t_onw + C_ind + 1
        ind2 <- t_inf + C_ind + 1
        ward_from <- data$ctd_timed_matrix[from, ind1]
        ward_to <- data$ctd_timed_matrix[to, ind2]
      }
      out <- tibble(from = from, to = to, ward_from = ward_from, ward_to = ward_to, kappa = kappa, t_inf = t_inf)
    }
    return(out)
  }

  ## summarise ward information across all MCMC steps
  get_summary <- function(i, res) {

    ## extract data
    nm <- names(res)
    from <- unlist(res[i, grep('alpha', nm)])
    to <- seq_along(from)
    t_inf <- unlist(res[i, grep('t_inf', nm)])
    t_onw <- unlist(res[i, grep('t_onw', nm)])
    kappa <- unlist(res[i, grep('kappa', nm)])

    if(is.null(t_onw)) t_onw <- rep(1, length(from))
    t_onw[is.na(t_onw)] <- -1

    ## calculate false negatives
    mapply(extract, from, to, t_inf, t_onw, kappa, SIMPLIFY = FALSE) %>%
      bind_rows()

  }

  ## assign some variables for later use
  data$ctd_timed_matrix <- data$ctd_timed_matrix[[1]]
  C_ind <- data$C_ind

  ## match ward names to incices used in outbreaker
  ward_match <- unique(data$ctd_timed[,2])[!is.na(unique(data$ctd_timed[,2]))]

  lapply(seq_len(nrow(res)), get_summary, res) %>%
    bind_rows() %>%
    mutate(
      across(c(ward_from, ward_to), ~ifelse(. == 0, "Unknown", ward_match[.])),
      across(c(from, to), ~data$ids[.])
    ) %>%
    arrange(to)

}

## this is a recursive function calculating various properties of the leaf node
##  -depth is the number of edges from a root node
##  -identify of the root node
## in the case of multiple parent nodes, the edge with the maximum value in
## 'rank_contact' is taken, which is either an edge attribute, or an edge
## attribute calculated from a node attribute by taking the difference in node
## attributes (e.g. difference in times of infection if rank_contact = 't_inf')
get_treestat <- function(i, depth, subtree_size, contacts, linelist, leaf,
                         rank_contact, reverse_rank_contact, leaf_parent) {

  ## identify parent by choosing highest ranked edge
  parents_ind <- which(contacts$to == i)
  parents <- contacts$from[parents_ind]
  if(length(parents) == 0) {
    parent <- NA
  } else if(length(parents) > 1) {
    edge_rank <- contacts[[rank_contact]][parents_ind]
    if(!reverse_rank_contact) {
      parent <- parents[which.max(edge_rank)]
    } else {
      parent <- parents[which.min(edge_rank)]
    }
  } else {
    parent <- parents
  }

  ## store parent of the leaf
  if(is.null(leaf_parent)) leaf_parent <- parent

  ## error if network is cyclical
  if(!is.na(parent) & parent == leaf) {
    stop("pruned network still contains cycles")
  }

  ## return results if root is hit or go one layer higher
  if(is.na(parent) || parent == 0) {
    return(list(depth = depth,
                subtree_size = subtree_size,
                parent = leaf_parent,
                root = i))
  } else {
    parent_ind <- linelist$id == parent
    get_treestat(i = parent,
                 depth = depth + 1,
                 replace(subtree_size, parent_ind, subtree_size[parent_ind] + 1),
                 contacts,
                 contacts,
                 leaf,
                 rank_contact,
                 reverse_rank_contact,
                 leaf_parent)
  }
}

## identify cycles in the network and remove the edge with the lowest rank in
## 'rank_contact'
prune_cycles <- function(i, leaf, contacts, cycle_edges,
                         rank_contact, reverse_rank_contact) {

  ## choose parent with highest ranked edge
  parents <- contacts[which(contacts$to == i),]
  if(nrow(parents) == 0) {
    return(contacts)
  } else if(nrow(parents) > 1) {
    if(!reverse_rank_contact) {
      parent <- parents[which.max(parents[[rank_contact]]),]
    } else {
      parent <- parents[which.min(parents[[rank_contact]]),]
    }
  } else {
    parent <- parents
  }

  ## return contacts if root found
  if(is.na(parent$from) || parent$from == 0) {
    return(contacts)
  }

  ## if this parent exists in the cycle, a cycle has been identified, the lowest
  ## ranked edge is removed and the pruning function is restarted to make sure no other
  ## cycles exist
  if(parent$from %in% cycle_edges$to) {
    cycle_edges <- rbind(cycle_edges, parent)
    edge_remove <- cycle_edges[which.min(cycle_edges[[rank_contact]]),]
    contacts <- contacts[-which(contacts$from == edge_remove$from &
                                contacts$to == edge_remove$to),]
    ## restart loop from leaf with pruned contacts
    contacts <- prune_cycles(leaf, leaf, contacts, NULL,
                             rank_contact, reverse_rank_contact)

  } else {

    ## if no loop, move onwards
    cycle_edges <- rbind(cycle_edges, parent)
    contacts <- prune_cycles(parent$from, leaf, contacts, cycle_edges,
                             rank_contact, reverse_rank_contact)

  }

  return(contacts)

}

## get_child_pos outputs children node positions relative to the parent
## position_dodge specifies if a parent and child can share the same y position
## parent_pos specifies where the parent sits relative to the children
get_child_pos <- function(n_children,
                          parent_pos = 'middle',
                          custom_parent_pos = NULL,
                          position_dodge = FALSE) {

  if(is.null(custom_parent_pos)) {

    if(parent_pos == 'middle') {

      if(n_children %% 2 != 0 & position_dodge) {
        n_children <- n_children + 1
        odd <- TRUE
      } else {
        odd <- FALSE
      }
      sq <- seq(1, 1 + 2*(n_children - 1), 2)
      out <- sq - stats::median(sq)
      if(odd & position_dodge) out <- out[-1]

    } else if(parent_pos == 'bottom') {

      if(position_dodge) {
        out <- 1:n_children
      } else {
        out <- 0:(n_children - 1)
      }

    } else if(parent_pos == 'top') {

      if(position_dodge) {
        out <- (1:n_children)*-1
      } else {
        out <- (0:(n_children - 1)*-1)
      }

    }
  } else {

    out <- custom_parent_pos(n_children)

  }

  return(out)

}

## this function will order the cases in such a manner that parent nodes are
## closest to their children, starting at the leaves of the tree and moving to
## the root. It returns the x and y coordinates of each node.
get_coor <- function(x,
                     x_axis,
                     position_dodge = FALSE,
                     root_order = 'subtree_size',
                     node_order = 'subtree_size',
                     reverse_root_order = FALSE,
                     reverse_node_order = FALSE,
                     rank_contact = x_axis,
                     reverse_rank_contact = FALSE,
                     unlinked_pos = 'bottom',
                     axis_type = c("single", "double", "none"),
                     parent_pos = 'middle',
                     custom_parent_pos = NULL,
                     method = 'ttree',
                     igraph_type = NULL) {

  axis_type <- match.arg(axis_type)

  ## add cluster membership, for use in root clustering
  x <- get_clusters(x)
  linelist <- x$linelist
  contacts <- x$contacts

  ## if rank_contact is a node attribute, calculate an edge attribute for each
  ## contact by taking the difference in node attributes
  if(!rank_contact %in% names(contacts)) {

    if(!rank_contact %in% names(linelist)) {
      stop("rank_contact is not found in linelist or contacts")
    }

    if(!inherits(linelist[[rank_contact]],
                 c("Date", "numeric", "integer", "POSIXct", "POSIXt"))) {
      stop("rank_contact must indicate a Date, numeric or integer value")
    }

    from_ind <- match(contacts$from, linelist$id)
    to_ind <- match(contacts$to, linelist$id)

    weight <- linelist[[rank_contact]][to_ind] - linelist[[rank_contact]][from_ind]
    contacts[[rank_contact]] <- as.numeric(weight)

  }

  N <- nrow(linelist)

  ## depth is the number of generations between a leaf and its root (+ 1)
  depth <- rep(0, N)

  ## parent is the parent node in the 'scaffold' tree (if multiple parents
  ## exist, one of these is selected to build the scaffold tree)
  parent <- root <- rep(NA, N)

  ## subtree size is the number of nodes downstream of a given node
  subtree_size <- rep(1, N)

  ## cycles are removed when building the scaffold tree by recursively removing
  ## the lowest ranked edge in a cycle (as defined by rank_contact)
  contacts_clean <- contacts
  for(i in linelist$id) {
    contacts_clean <- prune_cycles(i,
                                   leaf = i,
                                   contacts = contacts_clean,
                                   cycle_edges = NULL,
                                   rank_contact = rank_contact,
                                   reverse_rank_contact = reverse_rank_contact)
  }

  ## calculate various tree statistics using scaffold tree
  for(i in linelist$id) {
    treestat <- get_treestat(i,
                             depth = 1,
                             subtree_size = subtree_size,
                             contacts = contacts_clean,
                             linelist = linelist,
                             leaf = i,
                             rank_contact = rank_contact,
                             reverse_rank_contact = reverse_rank_contact,
                             leaf_parent = NULL)
    ind <- which(linelist$id == i)
    depth[ind] <- treestat$depth
    parent[ind] <- treestat$parent
    root[ind] <- treestat$root
    subtree_size <- treestat$subtree_size
  }

  ## add subtree size to linelist so that it can be called in node_order
  linelist$subtree_size <- subtree_size

  ## NAs treated as 0 to be safe, though these should be removed beforehand
  contacts$from[is.na(contacts$from)] <- 0

  ## pos is a matrix with the ranking of a node relative to its parent at each
  ## depth. every node inherits the ranking of its parent, ensuring nodes only
  ## split from their parents at the correct depth
  pos <- matrix(0, nrow(linelist), max(depth))

  ## iterate across every depth and split cases at that depth
  for(i in seq_len(max(depth))) {

    ## index of the cases found at a given depth i
    ind <- which(depth == i)

    ## child_pos_add will be the ranking of a child relative to its parent
    child_pos_add <- numeric(length(ind))

    ## non-root cases
    if(i > 1) {

      ## group nodes into siblings with shared parent
      siblings_grouped <- split(ind, parent[ind])

      for(siblings in siblings_grouped) {

        ## get the child position relative to the parent
        child_pos <- get_child_pos(length(siblings), parent_pos, custom_parent_pos, position_dodge)

        ## re-order siblings by order_nodes
        if(!is.null(node_order)) {
          siblings <- siblings[order(linelist[[node_order]][siblings],
                                     decreasing = reverse_node_order)]
        }

        ## add children position according to ordering
        child_pos_add[match(siblings, ind)] <- child_pos

      }

      ## children adopt ranking values from parent
      pos[ind, 1:(i-1)] <- pos[match(parent[ind], linelist$id), 1:(i-1)]

    } else {

      ## this section links roots that are connected by secondary connections
      n_parents <- table(contacts$to)
      multiple_parents <- names(n_parents)[n_parents > 1]
      if(length(multiple_parents) > 0) {

        ## identify cluster, group roots by cluster
        linelist$cluster_member <- as.numeric(linelist$cluster_member)
        root_ind <- match(unique(root), linelist$id)
        linked <- split(linelist$id[root_ind],
                        linelist$cluster_member[root_ind])

        ## get all combinations of linked roots
        get_combn <- function(j) if(length(unique(j)) > 1) t(utils::combn(unique(j), 2))
        linked_roots <- lapply(linked, get_combn)
        linked_roots <- linked_roots[!vapply(linked_roots, is.null, TRUE)]
        if(length(linked_roots) > 0) {
          linked_roots <- do.call(rbind, linked_roots)
        }

      } else {
        linked_roots <- numeric()
      }

      ## get position of roots
      child_pos <- get_child_pos(length(ind), parent_pos, custom_parent_pos, position_dodge)

      if(!is.null(root_order)) {

        ## get initial ordering of roots by root_order
        ord <- order(linelist[[root_order]][ind], decreasing = reverse_root_order)

        ## if linked roots exist, these are placed next to each other. a 'root
        ## of roots' which is highest ranked in root_order is chosen
        if(length(linked_roots) > 0) {

          ## bring the lower ranked root (from root_order) right below the
          ## higher ranked root by subtracting 0.0001
          for(j in seq_len(nrow(linked_roots))) {
            sub_ord <- ord[match(linked_roots[j,], linelist$id[ind])]
            min_ind <- match(linked_roots[j, which.min(sub_ord)], linelist$id[ind])
            max_ind <- match(linked_roots[j, which.max(sub_ord)], linelist$id[ind])
            ord[min_ind] <- ord[max_ind] - 0.0001*j
          }

          ## roots are now re-ordered the roots, considering only the 'roots of
          ## roots', ie the base root of a group of linked roots
          rr <- ord %in% ceiling(ord)

          ## this points to the root of each root
          rr_ind <- match(ceiling(ord), ord[rr])

          ## this gets the rank of the root of roots
          rr_ord <- rank(linelist[[root_order]][ind][rr],
                         ties.method = "first")

          ## this places the roots under the roots of roots
          diff <- ord - ceiling(ord)
          ord <- rr_ord[rr_ind] + diff

          ## get indices and get child positions
          ord <- match(ord, sort(ord))
          child_pos_add <- child_pos[ord]

        } else {

          ## add normal child positions if we don't have linked roots
          child_pos_add[ord] <- child_pos

        }

      } else {

        ## if node order is null, use native ordering in linelist
        child_pos_add <- child_pos

      }

    }

    ## add the child positions at the given depth
    pos[ind, i] <- child_pos_add

  }

  ## this is true if i is ranked higher than j (essentially looks at the
  ## difference in ranking at each depth and choose the first difference)
  is_higher <- function(i, j) {
    comp <- pos[i,] - pos[j,]
    out <- comp[comp != 0][1] > 0
  }

  ## rank all cases by pairwise comparison - this is quite slow for large
  ## networks, there is definitely a quicker way
  comb <- expand.grid(seq_len(nrow(pos)), seq_len(nrow(pos)))
  rel_pos <- matrix(FALSE, nrow = nrow(pos), ncol = nrow(pos))

  rel_pos[as.matrix(comb)] <- apply(comb, 1, function(x) is_higher(x[1], x[2]))
  rel_pos[is.na(rel_pos)] <- FALSE

  ## get global ranking by taking the sum of all pairwise rankings
  abs_pos <- apply(rel_pos, 1, sum)

  ## Get isolated cases and place them as specified by unlinked_pos
  ## If position_dodge, order these by root_order
  contacts <- contacts[contacts$from != 0,]
  unlinked <- which(!linelist$id %in% c(contacts$from, contacts$to))
  if(length(unlinked) > 0) {
    if(position_dodge) {
      unlinked_order <- order(linelist[[root_order]][unlinked],
                              decreasing = reverse_root_order)
      unlinked_to_add <- seq_along(unlinked)
    } else {
      unlinked_order <- seq_along(unlinked)
      unlinked_to_add <- 0
    }

    ## Place nodes at top or bottom by adding 10000 to abs_pos
    if(unlinked_pos == 'top') {
      abs_pos[unlinked[unlinked_order]] <- 10000 + unlinked_to_add
    } else if(unlinked_pos == 'bottom') {
      abs_pos[unlinked[unlinked_order]] <- -10000 + unlinked_to_add
    }
  }

  ## use igraph coordinate calculations instead
  if(!is.null(igraph_type)) {

    net <- igraph::graph_from_data_frame(stats::na.omit(x$contacts),
                                         vertices = x$linelist)

    if(igraph_type == 'rt') {
      abs_pos <- igraph::layout.reingold.tilford(net, root = which(depth == 1))[,1]
    } else if(igraph_type == 'sugiyama') {
      abs_pos <- igraph::layout.sugiyama(net, layers = x$linelist[[x_axis]])$layout[,1]
    } else if(igraph_type == 'fr') {
      abs_pos <- igraph::layout.fruchterman.reingold(net,
                                                 minx = x$linelist[[x_axis]],
                                                 maxx = x$linelist[[x_axis]])[,2]
    } else {
      stop("igraph_type must be one of 'rt', 'sugiyama' or 'fr'")
    }

  }

  ## this orders the cases and adds y_adj positions at the top and bottom to
  ## provide space for the axes
  y_adj <- 2
  if(method == 'ttree') {
    if(axis_type == 'double') {
      y_pos <- match(abs_pos, sort(unique(abs_pos)))
      y_pos = rescale(c(min(y_pos) - y_adj, max(y_pos) + y_adj, y_pos), 0, 1)
    } else if(axis_type == 'single') {
      y_pos <- match(abs_pos, sort(unique(abs_pos)))
      y_pos = rescale(c(min(y_pos) - y_adj, y_pos), 0, 1)
    } else if(axis_type == 'none') {
      y_pos <- match(abs_pos, sort(unique(abs_pos)))
      y_pos = rescale(y_pos, 0, 1)
    }
  } else {
    y_pos <- match(abs_pos, sort(unique(abs_pos)))
    y_pos <- rescale(y_pos, 0, 1)
  }

  ## Also return parent because we need scaffold tree for later
  return(list(y = y_pos, subtree_size = subtree_size))

}

## Rescale a vector of numerics to min and max
rescale <- function(x, min_val = 0, max_val = 1) {
  if(length(unique(x)) == 1L) return(x)
  out <- (x - min(x, na.rm = TRUE))/(diff(range(x, na.rm = TRUE)))*(max_val - min_val)
  out <- out + min_val
  return(out)
}

## Minimum genetic distance to another case for each case
get_mindist <- function(dist) {
  sapply(seq_len(nrow(dist)), function(x) min(dist[x, -x]))
}

## Fill in an MCMC state from a given row in the outbreaker output
update_par <- function(param, i, res, data, config) {

  unl <- function(...) unlist(..., use.names = FALSE)
  param$alpha <- as.vector(unl(res[i, grep('alpha', names(res))]))
  param$t_inf <- as.vector(unl(res[i, grep('t_inf', names(res))]))
  param$t_onw <- as.vector(unl(res[i, grep('t_onw', names(res))]))
  param$kappa <- as.integer(as.vector(unl(res[i, grep('kappa', names(res))])))

  param$pi <- unl(res[i, grep('pi', names(res))])
  param$eps <- unl(res[i, grep('eps', names(res))])
  param$lambda <- unl(res[i, grep('lambda', names(res))])
  param$eta <- unl(res[i, grep('eta', names(res))])
  param$tau <- unl(res[i, grep('tau', names(res))])
  param$mu <- unl(res[i, grep('mu', names(res))])

  param$t_onw[is.na(param$t_onw)] <- -1000
  param$kappa[is.na(param$kappa)] <- 1

  param <- update_trans(data, param, config)

  return(param)

}

## Update transition matrices
update_trans <- function(data, param, config) {
  for(i in seq_along(param$trans_mat)) {
    param$trans_mat[[i]] <- outbreaker2:::get_transition_mat(
      data$pp_trans_adj[[i]],
      data$pp_place[[i]],
      data$pp_place_adj[[i]],
      param$eps[i],
      param$tau[i],
      data$N_place_unobserved[i],
      config$max_kappa
    )
  }
  return(param)
}

## get minimum genetic distance across all posterior draws
get_res_mindist <- function(res, dist) {
  sapply(
    seq_along(res[grep("alpha", names(res))]),
    function(i) min(dist[i, unique(res[grep("alpha", names(res))][[i]])])
  )
}

## get potential ancestors for case i
get_potential_old <- function(i, dist, dates) {
  mindist <- get_mindist(dist)[i]
  potential <- which(dist[i,] == mindist)
  potential <- potential[potential != i]
  ## get those with onsets before i
  potential <- potential[data$dates[i] + 5 > data$dates[potential]]
  if(length(potential) == 0) return(numeric())
  potential[sapply(
    potential,
    function(x)
      any(get_ward(x, (data$dates[i]-30):(data$dates[i] - 1), data) != -1)
  )]
}

## extract transmission pair information
get_tpairs <- function(res, data) {
  lapply(
    seq_len(sum(grepl("alpha", names(res)))),
    function(i) {
      tibble(
        from = res[[glue("alpha_{i}")]],
        to = i,
        t_inf = res[[glue("t_inf_{i}")]],
        t_onw = res[[glue("t_onw_{i}")]],
        kappa = res[[glue("kappa_{i}")]],
        ward_from = mapply(get_ward, from, ifelse(is.na(kappa) | kappa == 1, t_inf, t_onw), MoreArgs = list(data = data)),
        ward_to = mapply(get_ward, to, t_inf, MoreArgs = list(data = data)),
        dist = data$D[cbind(from, to)]
      )
    }
  ) %>%
    bind_rows() %>%
    drop_na(from)
}

## get potential ancestors on the same ward with a max genetic distance
get_potential <- function(i, data, max_dist = 0) {

  onset <- data$dates[i]
  dates <- (onset - 30):(onset - 1)
  potential <- which(data$D[i,] <= max_dist)
  potential <- potential[potential != i]

  wards <- get_ward(i, dates, data)
  potential_wards <- sapply(potential, get_ward, dates, data)

  ## find the case-date where the ward is the same
  matches <- apply(potential_wards, 2, function(x) x == wards & !x %in% c(-1, 0))
  ind <- which(matches, arr.ind = TRUE)

  tibble(
    from = potential[ind[,2]],
    to = i,
    t_inf = dates[ind[,1]],
    t_ons_to = data$dates[i],
    t_ons_from = data$dates[from],
    serial_interval = t_ons_to - t_ons_from,
    incub_to = t_ons_to - t_inf,
    ward_to = mapply(get_ward, to, t_inf, MoreArgs = list(data = data)),
    ward_from = mapply(get_ward, from, t_inf, MoreArgs = list(data = data))
  ) %>%
    filter(t_inf < t_ons_from)

}

## get the ward likelihood from the wards
get_ward_ll <- function(w1, w2, trans_mat, kappa = 1) {
  arr <- array(param$trans_mat[[1]], dim = c(100, 100, 5))
  arr[w1 + 1, w2 + 1, kappa]
}

## poisson ll with i subsetting
ll_pois_i <- function(data, param, i = NULL) {
  if(is.null(i)) i <- seq_along(param$alpha)
  sum(
    dpois(
      data$D[cbind(param$alpha[i], i)],
      mu*param$kappa[i],
      log = TRUE
    ),
    na.rm = TRUE
  )
}

## Get alpha and beta of beta distribution from mean and var
get_beta_par <- function(mu, sd = NULL, cv = NULL) {
  if(!is.null(sd)) var <- sd^2
  if(is.null(var)) var <- (mu*cv)^2
  if(var >= mu*(1-mu)) stop("var is too large")
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(c(alpha, beta))
}

## check if there is any overlap between wards stays
get_overlap <- function(start, end) {
  map_lgl(
    seq_along(start),
    function(i) any(start[i] > end[-i] & start[-i] > end[i])
  )
}

## scale dates up (i.e. expand) or down (i.e. contract) by a factor fiven by scale
scale_ward_dates <- function(ward, origin, to = c("up", "down"), scale = 6) {
  to <- match.arg(to)
  if(to == "up") {
    mutate(
      ward,
      across(
        c(start, end),
        ~ (as.numeric(.x) - as.numeric(origin)) %>%
          `*`(scale) %>%
          as.POSIXct(origin = origin, tz = "GMT") %>%
          floor_date("day") %>%
          as.Date()
      )
    ) %>%
      group_by(id, start) %>%
      ## this will shuffle dates by +/- days so we capture all wards
      do(correct_dates(.)) %>%
      ## as this works recursively we have do it multiple times - twice seems to
      ## be enough to ensure dates are shuffled enough to capture all ward stays
      do(correct_dates(.)) %>%
      ungroup()
  } else {
    mutate(
      ward,
      across(
        c(start, end),
        ~ as.POSIXct(with_tz(.x, "GMT")) %>%
          {as.numeric(.) - as.numeric(origin)} %>%
          `/`(scale) %>%
          as.POSIXct(origin = origin, tz = "GMT")
      )
    )
  }
}

## scale dates up (i.e. expand) or down (i.e. contract) by a factor fiven by scale
scale_onset_dates <- function(onset, origin, to = c("up", "down"), scale = 6) {
  to <- match.arg(to)
  if(to == "up") {
    as.POSIXct(with_tz(onset, "GMT"), format = "%Y-%m-%d %H:%M:%S", tz = "GMT") %>%
    {as.numeric(.) - as.numeric(origin)} %>%
      `*`(scale) %>%
      as.POSIXct(origin = origin, tz = "GMT") %>%
      floor_date("day") %>%
      as.Date()
  } else {
    as.POSIXct(with_tz(onset, "GMT")) %>%
      {as.numeric(.) - as.numeric(origin)} %>%
      `/`(scale) %>%
      as.POSIXct(origin = origin, tz = "GMT")
  }
}

## shuffle dates so that all ward event can be captured
correct_dates <- function(df) {
  if(nrow(df) > 1) {
    df %<>% arrange(end)
    ## this is the event happening only on that "day"
    single <- df$start - df$end == 0
    if(sum(single) > 1)
      warning(paste0("resolution is too small to include all data points; ",
              "if you see many of these warnings, increase the scale"))
    for(i in which(!single)) {
      ## if before
      if(i < which(single)[1]) df$end[i] <- df$end[i] - 1
      if(i > which(single)[1]) df$start[i] <- df$start[i] + 1
    }
  }
  return(df)
}

## rescale distributions
scale_distribution <- function(density, scale = 6) rep(density/scale, each = scale)

## calculate the individual likelihood for all cases for a given step in the MCMC
get_step_likelihood <- function(step, res, data, param, config, likelihoods) {
  param <- update_par(param, step, res, data, config)
  ll <- map_dbl(seq_along(data$dates), ~ outbreaker2:::cpp_ll_all(data, param, .x, likelihoods))
  ## we don't calculate likelihoods for imports as they are only calculated on
  ## onset dates and aren't comparable with other likelihoods
  ll[is.na(param$alpha)] <- NA
  return(ll)
}

## calculate the average likelihood across the entire run
get_average_likelihood <- function(res, data, param, config, likelihoods, burnin = 1000) {

map_dfc(
    which(res$step > burnin),
    get_step_likelihood,
    res, data, param, config, likelihoods
  ) %>%
    apply(1, mean, na.rm = TRUE)

}

## assign imports according the import probability and the case
## likelihoods. n_chains specifies how many different outbreaker chains you want
## to run with different imports draws.
assign_imports <- function(n_chains, p_import_mu, p_import_sd, res, ll, config) {

  ## these are fixed as imports from the previous run
  prop_fixed <- mean(is.na(config$init_tree))

  ## check that p_import_mu is high enough
  if(p_import_mu < prop_fixed) {
    stop(paste0("p_import_mu is lower than the minimum necessary number of imports",
                " given the ward data and onset dates"))
  }

  ## get beta prior values for our import probability - we need to calculate the
  ## probability for the unfixed imports in a manner that ensure unfixed + fixed
  ## = p_import_mu
  par <- get_beta_par((p_import_mu - prop_fixed)/(1 - prop_fixed), p_import_sd)

  ## draw import probabilities
  p_import <- rbeta(n_chains, par[1], par[2])

  ## draw non-fixed import numbers
  n_import <- map_dbl(p_import, ~ rbinom(1, length(ll) - sum(is.na(config$init_tree)), .x))

  ## cases that are not fixed as imports
  unfixed <- !is.na(config$init_tree)

  ## draw import indices from unfixed (with 1/ll as the probability of being an import)
  imports <- map(n_import, ~ sample(which(unfixed), .x, FALSE, exp(-ll[unfixed])))

  ## set last tree of res as first tree in config
  param <- update_par(create_param()$current, nrow(res), res, data, config)

  ## assign
  config$init_tree <- config$init_alpha <- param$alpha
  config$init_t_inf <- param$t_inf
  config$init_t_onw <- param$t_onw
  config$init_kappa <- param$kappa
  config$init_pi <- param$pi
  config$init_mu <- param$mu
  config$init_eps <- param$eps
  config$init_tau <- param$tau

  ## assign these indices as imports and return
  map(imports, ~ modify_in(config, "init_tree", modify_at, .x, ~ NA))

}

## run outbreaker with a given config and save as run_i
run_and_save <- function(config, i, data, priors, likelihoods) {
  res <- outbreaker(data, config, priors = priors, likelihoods = likelihoods)
  export(res, here(glue("outputs/run_{i}.rds")))
  return(res)
}
