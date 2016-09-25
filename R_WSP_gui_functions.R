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
  assign('exp_set', es_tmp, .GlobalEnv);
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
      #testing only
      em = exprs(exp_set);
      half = floor(dim(em)[2]);
      em1 = em[,1:half];
      em2 = em[,-(1:half)];
      diff = gene_diff(em1, em2, svalue(pval_slider), svalue(FC_slider), svalue(outdir_edit));
      assign('f_res', diff, .GlobalEnv);
    }else if(user_choice == "Clustering"){
      Klasteryzacja(exp_set, svalue(clnum_slider), svalue(outdir_edit));
    }else if(user_choice == "Clustering 2"){
      clusters = klast(exp_set, svalue(outdir_edit));
      assing('f_res', clusters, .GlobalEnv);
    }
  }
  addHandlerClicked(run_button, handler = run_analysis);
  return();
}