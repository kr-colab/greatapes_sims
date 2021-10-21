build_nwk = function(focal, t, parent_col="parent", child_col="edge", len_col= "gens") {
    #print(paste0("got in: ", focal))
    treestr = paste0(paste0(focal,":",t[t[child_col]==focal,][len_col]))
    if (t[t[child_col]==focal,][parent_col]=="") {
        treestr = paste0(treestr, ";")
    }
    #print(new_bit)
    childs = t[t[parent_col]==focal,]
    if (nrow(childs) == 0) {
        return (treestr)
    }
    else {
        stopifnot(nrow(childs) == 2)
        treestr = paste0("(",build_nwk(childs[child_col][1,], t),",", build_nwk(childs[child_col][2,], t),")", treestr)
    }
    return (treestr)
}

harmonic_ne = function(focal) {
    ## THIS DOESN"T MAKE MUCH SENSE?
    path = nodepath(tree, treetbl[treetbl$label==root,]$node,treetbl[treetbl$label==focal,]$node)
    treetbl[treetbl$node %in% path,]$label
    times = treetbl[treetbl$node %in% path,]$branch.length
    t = sum(treetbl[treetbl$node %in% path,]$branch.length)

    ns = treetbl[treetbl$node %in% path,]$N
    inv_ne = sum((times/t)*(1/ns))
    return(1/inv_ne)
}

id_to_label = function(nid, treetbl) {
    return(treetbl[treetbl$node == nid,]$label)
}

label_to_id = function(lab, treetbl) {
    return(treetbl[treetbl$label == lab,]$node)
}

path_between_tips = function(tip1, tip2, mrca, treetbl, tree) {
    # FINDS PATH BETWEEN TIPS IN THE TREE
    path = sapply(nodepath(tree, label_to_id(tip1, treetbl), label_to_id(tip2, treetbl)), id_to_label, treetbl=treetbl)
    path = setdiff(path, mrca)
    external_branches = tree$tip.label
    return(list(external = path[path %in% external_branches], internal = path[!path %in% external_branches]))
}

get_mrow = function(row, cols, tree, treetbl) {
    # GETS THE COEFFICIENTS FOR THE NEW VARS (COLS) FROM THE OLD VAR IN ROW
    mrow = rep(0, length(cols))
    path = path_between_tips(row[1], row[2], row[4], treetbl, tree)
    if(row[7] == "pi") {
        mrow[cols == paste("pi", row[1], sep="_")] = 1
    } else {
        # dxys consist of the tree height in the tips (pi/2) + the sum of
        # branch components (new variables called N)
        #print(path)
        mrow[cols %in% paste("pi", path$external, sep="_")] = 1/2
        mrow[cols %in% paste("N", Reduce(c, path), sep="_")] = 1
    }
    return(mrow)
}

meta_from_fname = function(fname, prop=NULL) {
    ga_data_str = "greatapes-diversity-data"
    is_ga_data = grepl(ga_data_str,inpath, fixed=TRUE)
    if (is_ga_data) {
        strp = '.+win-size_(\\d+)_merged-mask_(\\w+)'
        if (is.null(prop)){
            strp = paste0(strp, '_prop-acc_(.+)')
        }
        strp = paste0(strp, "\\.tsv")
        matches = str_match(fname, strp)
        win_size = matches[2]
        merged_mask = matches[3]
        if (is.null(prop)) {
            prop = matches[4]
        }
        spaced_desc = paste0("win-size=", win_size, " merged-mask=", merged_mask, " prop-acc=", prop)
        desc = str_replace_all(spaced_desc, " ", "_")
        desc = str_replace_all(desc, "=", "_")
        meta = list("win_size" = as.integer(win_size), "merged_mask" = merged_mask, "spaced_desc"=spaced_desc, "desc" = desc, "prop" = as.numeric(prop), "is_ga_data"=is_ga_data)
    } else {
        strp = '.+sup-rand-id_(.+)_rep_(\\d+)_win-size_(\\d+)_sample-size_(\\d+)\\.tsv'
        matches = str_match(fname, strp)
        suprand = matches[2]
        rep = matches[3]
        win_size = matches[4]
        sample_size = matches[5]
        spaced_desc = paste0("sup-rand-id=", suprand, " rep=", rep, "\nwin-size=", win_size, " sample-size=", sample_size, " prop-acc=", prop)
        desc = str_replace_all(spaced_desc, " ", "_")
        desc = str_replace_all(desc, "=", "_")
        desc = str_replace_all(desc, "\n", "_")
        meta = list("win_size" = as.integer(win_size), "sup_rand_id"= suprand, "rep"=as.integer(rep), "sample_size"=as.integer(sample_size), "spaced_desc"=spaced_desc, "desc" = desc, "prop" = as.numeric(prop), "is_ga_data"=is_ga_data)
    }
    return (meta)
}
