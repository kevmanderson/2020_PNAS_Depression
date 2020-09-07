library(tidyverse)
library(gifti)
library(XML)


map2color = function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}


plot_label = function(plot_vals, limits, mypal, parcel_num, net_num, out_path){

    base_dir         = '/gpfs/milgram/project/holmes/kma52/mdd_sst'
    dlabel_template  = paste0(base_dir, '/reference_files/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order.dscalar.nii')
    dscalar_template = paste0(base_dir, '/reference_files/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order.dlabel.nii')
    parcel_info_file = paste0(base_dir, '/reference_files/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order_info.txt')

    label_in = read_csv(parcel_info_file, col_names=F)
    label_array = label_in$X1[grep('Net', label_in$X1)]

    rgb_matrix = t(col2rgb(map2color(plot_vals, mypal, limits=c(limits[1], limits[2]) )))/256
    rgb_matrix[is.na(plot_vals),] =  182/256

    # write label color reference file
    write_text_file = gsub('.dlabel.nii', '.txt', plot_dlabel)
    write_labels    = NULL
    for (row in 1:nrow(rgb_matrix)){
        label        = label_array[row]
        write_rgb    = paste(rgb_matrix[row,]*256, collapse=' ')
        write_row    = paste(as.character(row), write_rgb, '1',  collapse=' ')
        write_labels = c(write_labels, label)
        write_labels = c(write_labels, write_row)
    }
    write.table(x=write_labels, file=write_text_file, col.names=F, row.names=F, quote=F)

    import_label_cmd = c('/gpfs/milgram/apps/hpc.rhel7/software/ConnectomeWorkbench/1.3.2/bin_rh_linux64/wb_command -cifti-label-import',dscalar_template, write_text_file, plot_dlabel)
    system(paste(import_label_cmd, collapse=' '))
}


# function to plot parcel-wise cortical correlations for viewing with HCP wb_view
plot_matlab = function(values, out_path, parcel_num, net_num){
    base_dir = '/gpfs/milgram/project/holmes/kma52/mdd_gene_expr'
    dscalar_template = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order.dscalar.nii')
    parcel_info_file = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order_info.txt')
    write_val_file   = paste0('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/ahba/tmp/tmp.txt')
    write_delim(delim=' ', x=as.data.frame(values), path=write_val_file,col_names=F)
    save_path = out_path
    print(save_path)

    matfunc = 'plotVolOnSurface'
    cmd     = paste0('/gpfs/milgram/apps/hpc.rhel7/software/MATLAB/2017b/bin/matlab -nosplash -nodesktop -r "cd(\'/gpfs/milgram/project/holmes/kma52/mdd_gene_expr/scripts/util\');',
                    matfunc, '(\'', dscalar_template, '\',\'', parcel_info_file, '\',\'',
                    write_val_file, '\',\'', save_path, '\'); exit;"')
    system(cmd)
}



plot_foci = function(mni.coords, vertices, foci.names, ref.surface, rgba.color, wb.path, projection.structure, out.name){

    gifti.file  <- readgii(ref.surface)


  top     <- newXMLNode("FociFile", attrs=c(Version="1.0"))
  metadat <- newXMLNode("MetaData", parent = top)
  md      <- newXMLNode("MD", parent = metadat)
  name    <- newXMLNode("Name", parent = md)
  value   <- newXMLNode("Name", parent = md)
  newXMLCDataNode('Caret-Version', parent = name, doc = NULL, at = NA, sep = "\n")
  newXMLCDataNode('5.64', parent = value, doc = NULL, at = NA, sep = "\n")

  md2      <- newXMLNode("MD", parent = metadat, sibling=md)
  name      <- newXMLNode("Name", parent = md2)
  value      <- newXMLNode("Value", parent = md2)
  newXMLCDataNode('Date', parent = name, doc = NULL, at = NA, sep = "\n")
  newXMLCDataNode('2017-05-14T16:03:20', parent = value, doc = NULL, at = NA, sep = "\n")

  md3      <- newXMLNode("MD", parent = metadat, sibling=md)
  name      <- newXMLNode("Name", parent = md3)
  value      <- newXMLNode("Value", parent = md3)
  newXMLCDataNode('UserName', parent = name, doc = NULL, at = NA, sep = "\n")
  newXMLCDataNode('vanessen', parent = value, doc = NULL, at = NA, sep = "\n")

  md4      <- newXMLNode("MD", parent = metadat, sibling=md)
  name      <- newXMLNode("Name", parent = md4)
  value      <- newXMLNode("Value", parent = md4)
  newXMLCDataNode('comment', parent = name, doc = NULL, at = NA, sep = "\n")
  newXMLCDataNode('none', parent = value, doc = NULL, at = NA, sep = "\n")

  labels <- newXMLNode("LabelTable", parent = top)
  #base.node <- newXMLNode("Label", parent=labels)
  node.list <- list()
  for (idx in 1:dim(rgba.color)[1]){

    # Hacky fix for a single vertex that gives us trouble
    vertex   <- vertices[idx]
    v.idxs   <- which(gifti.file$data$triangle == vertex, arr.ind=TRUE)
    tri.vertices <- gifti.file$data$triangle[v.idxs[,1],]
    tmp.vert     <- tri.vertices[1:2,]
    vert.six <- c(tmp.vert[1,], tmp.vert[2,])
    TriAnat <- NULL
    for (i in 1:6){
      tmp <- gifti.file$data$pointset[vert.six[i],]
      TriAnat <- c(TriAnat, tmp)
    }
    if (length(TriAnat) != 18){
      write(idx,'')
      next
    }

    cur.rgba  <- rgba.color[idx,]
    cur.label <- paste('Label_', foci.names[idx], sep='')
    tmp.node <- newXMLNode("Label", attrs=c(Key=idx, Red=cur.rgba$r, Green=cur.rgba$g, Blue=cur.rgba$b, Alpha=cur.rgba$a))
    newXMLCDataNode(cur.label, parent=tmp.node, doc = NULL, at = NA, sep = "\n")
    node.list[[idx]] <- tmp.node
  }
  addChildren(labels, node.list)

  foci.list <- list()
  for (idx in 1:dim(rgba.color)[1]){
    cur.label <- paste('Label_', foci.names[idx], sep='')
    cur.mni   <- mni.coords[idx,]
    foci.node <- newXMLNode('Focus', attrs=c(Index=idx))
    name.node <- newXMLNode('Name')
    newXMLCDataNode(cur.label, parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    xyz.node <- newXMLNode('SearchXYZ',parent=foci.node)
    newXMLTextNode(paste(cur.mni$mni_x, cur.mni$mni_y, cur.mni$mni_z), parent=xyz.node)
    name.node <- newXMLNode('Geography')
    newXMLCDataNode('', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('Area')
    newXMLCDataNode('', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('RegionOfInterest')
    newXMLCDataNode('', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('Size')
    newXMLCDataNode('0', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('Statistic')
    newXMLCDataNode('', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('Comment')
    newXMLCDataNode('', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('ClassName')
    newXMLCDataNode('', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('SumsIDNumber')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('SumsRepeatNumber')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('SumsParentFocusBaseID')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('SumsVersionNumber')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('SumsMSLID')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('AttributeID')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('StudyMetaDataLinkSet')
    name.node <- newXMLNode('AttributeID')
    newXMLCDataNode('-1', parent=name.node, doc = NULL, at = NA, sep = "\n")
    addChildren(foci.node, name.node)
    name.node <- newXMLNode('SurfaceProjectedItem')
    struct.node <- newXMLNode('Structure', parent=name.node)
    newXMLTextNode(projection.structure, parent=struct.node)
    struct.node <- newXMLNode('StereotaxicXYZ', parent=name.node)
    newXMLTextNode(paste(cur.mni$mni_x, cur.mni$mni_y, cur.mni$mni_z), parent=struct.node)
    struct.node <- newXMLNode('VolumeXYZ', parent=name.node)
    newXMLTextNode(paste(cur.mni$mni_x, cur.mni$mni_y, cur.mni$mni_z), parent=struct.node)
    addChildren(foci.node, name.node)

    vertex   <- vertices[idx]
    v.idxs   <- which(gifti.file$data$triangle == vertex, arr.ind=TRUE)
    tri.vertices <- gifti.file$data$triangle[v.idxs[,1],]
    tmp.vert     <- tri.vertices[1:2,]
    vert.six <- c(tmp.vert[1,], tmp.vert[2,])
    TriAnat <- NULL
    for (i in 1:6){
      tmp <- gifti.file$data$pointset[vert.six[i],]
      TriAnat <- c(TriAnat, tmp)
    }
    if (length(TriAnat) != 18){
      write(idx,'')
      next
    }
    ve.node <- newXMLNode('VanEssenProjection')

    struct.node <- newXMLNode('DR', parent=ve.node)
    newXMLTextNode(0, parent=struct.node)
    struct.node <- newXMLNode('TriAnatomical', parent=ve.node)
    newXMLTextNode(paste(TriAnat, collapse=' '), parent=struct.node)
    struct.node <- newXMLNode('ThetaR', parent=ve.node)
    newXMLTextNode(0, parent=struct.node)
    struct.node <- newXMLNode('PhiR', parent=ve.node)
    newXMLTextNode(0, parent=struct.node)
    struct.node <- newXMLNode('TriVertices', parent=ve.node)
    newXMLTextNode(paste(vert.six, collapse=' '), parent=struct.node)
    two.vert <- names(rev(sort(table(vert.six))))[1:2]
    struct.node <- newXMLNode('Vertex', parent=ve.node)
    newXMLTextNode(paste(two.vert, collapse=' '), parent=struct.node)
    struct.node <- newXMLNode('VertexAnatomical', parent=ve.node)
    tmp1 <- gifti.file$data$pointset[as.numeric(two.vert[1]),]
    tmp2 <- gifti.file$data$pointset[as.numeric(two.vert[2]),]
    newXMLTextNode(paste(c(tmp2, tmp1), collapse=' '), parent=struct.node)
    struct.node <- newXMLNode('PosAnatomical', parent=ve.node)
    newXMLTextNode(paste(cur.mni$mni_x, cur.mni$mni_y, cur.mni$mni_z), parent=struct.node)
    struct.node <- newXMLNode('FracRI', parent=ve.node)
    newXMLTextNode(5, parent=struct.node)
    struct.node <- newXMLNode('FracRJ', parent=ve.node)
    newXMLTextNode(5, parent=struct.node)

    addChildren(name.node, ve.node)
    addChildren(foci.node, name.node)

    foci.list[[idx]] <- foci.node
  }
  addChildren(top, foci.list)
  saveXML(top, file=out.name)
}





