get_file <- function(h, ...){
  dir = file.choose();
  svalue(h$action) = dir;
  assign("data_file_dir", dir, .GlobalEnv);
  return();
}

preprocess_CEL_data <- function(h, ...){
  my_wd = getwd();
  dfd = get('data_file_dir', envir = .GlobalEnv);
  setwd(dirname(dfd));
  es_tmp = get_expset(dfd);
  #es_cust = klast(es_tmp);removed
  assign('exp_set', es_tmp, .GlobalEnv);
  #assign('clusters', es_clust, .GlobalEnv);removed
  setwd(my_wd);
  return();
}

create_tools <- function(h, ...){
  user_choice = svalue(analysis_radio);
  #recreate boxes
  if(exists("an_settings", envir = .GlobalEnv)){
    delete(settings_container, an_settings);
  }
  an_settings = gframe(text = "Analysis Settings", horizontal = F, container = settings_container);
  assign('an_settings', an_settings, .GlobalEnv);
  if(exists("an_output_settings", envir = .GlobalEnv)){
    delete(settings_container, an_output_settings);
  }
  an_output_settings = gframe(text = "Output Settings", horizontal = F, container = settings_container);
  assign('an_output_settings', an_output_settings, .GlobalEnv);
  if(exists('settings_list', envir = .GlobalEnv)){
    rm(settings_list, pos = globalenv());
  }
  #fill boxes with tools
  if(user_choice == "Genes Difference"){
    glabel(text = "FC", container = an_settings);
    pval_slider = gslider(from = 0.01, to = 0.1, by = 0.01, value = 0.05, container = an_settings);
    glabel(text = "p-val threshold", container = an_settings);
    FC_slider = gslider(from = -7, to = 7, by = 0.5, value = 0, container = an_settings);
    setlist = list(pval_slider, FC_slider);
    names(setlist) = c('pval', 'fc');
    assign('settings_list', setlist, envir = .GlobalEnv);
  }else if(user_choice == "Clustering"){
    glabel("Number of clusters", container = an_settings);
    clnum_slider = gslider(from = 2, to = 10, by = 1, value = 3, container = an_settings);
    setlist = list(clnum_slider);
    names(setlist) = c('clustnum');
    assign('settings_list', setlist, envir = .GlobalEnv);
  }
  #create out control
  outdir_edit = gedit("Results save directory", container = an_output_settings);
  outdir_choose = gbutton('Choose directory', container = an_output_settings);
  out_choose <- function(h, ...){
    dir = choose.dir(default = getwd());
    svalue(h$action) = dir;
    return();
  }
  addHandlerClicked(outdir_choose, handler = out_choose, action = outdir_edit);
  #create run control
  run_button = gbutton("Run analysis!", container = an_output_settings);
  #create run function
  run_analysis <- function(h, ...){
    svalue(outdir_edit) = gsub("\\\\","/", svalue(outdir_edit));
    exp_set = get('exp_set', envir = .GlobalEnv);
    if(user_choice == "Genes Difference"){
      #run this
    }else if(user_choice == "Clustering"){
      Klasteryzacja(exp_set, svalue(clnum_slider));
    }else if(user_choice == "Clustering 2"){
      clusters = klast(exp_set, svalue(outdir_edit));
    }
  }
  addHandlerClicked(run_button, handler = run_analysis);
  return();
}

create_tools_gui <- function(){#obsolete?
  #gene diff
  gene_diff_group_input = ggroup(horizontal = F);
  file1_ged = gedit("File 1", container = gene_diff_group_input);
  file2_ged = gedit("File 2", container = gene_diff_group_input);
  pval_slider = gslider(from = 0.01, to = 0.05, by = 0.01, value = 0.05, container = gene_diff_group_input);
  FC_slider = gslider(from = 0.01, to = 0.1, by = 0.01, value = 0.05, container = gene_diff_group_input);
  gene_diff_group_output = ggroup(horizontal = F);
  gene_diff_in = list(gene_diff_group_input, file1_ged, file2_ged, pval_slider, FC_slider);
  gene_diff_out = list();
  #pca
  pca_group_input = ggroup(horizontal = F);
  file_ged = gedit("File 1", container = gene_diff_group_input);
  clnum_ged = gedit("Clusters", container = gene_diff_group_input);
  pca_group_output = ggroup(horizontal = T);
  pca_win1 = ggraphics(container = pca_group_output);
  pca_win2 = ggraphics(container = pca_group_output);
  pca_win3 = ggraphics(container = pca_group_output);
  pca_in = list(pca_group_input, file_ged, clnum_ged);
  pca_out = list(pca_group_output, pca_win1, pca_win2, pca_win3);
  elements_list = list(gene_diff_in, gene_diff_out, pca_in, pca_out);
  names(elements_list) = c("gene_diff_in", "gene_diff_out", "pca_in", "pca_out");
  return(elements_list); 
}