#>>>>>>>
## Example code from Gergo Dioszegi (Untested yet)
## Invoking CloudCompare command lines from R

#> Packages
require(lidR);require(tidyverse)

#> 1)
CC <- function(cc_function) {
CC_cmd = str_c(CC_dir, cc_function, sep = " ")
for (f in 1:length(CC_cmd)) {
system(command = CC_cmd[f])
}
}

#> 2)
cc_TREEiso <- function(file,
K1 = 5, L1 = 1.0, DEC_R1 = .05, # Cut-pursuit stage one
K2 = 20, L2 = 20, MAX_GAP = 2, DEC_R2 = 0.1, # Cut-pursuit stage two
VER_O_W = .5, RHO = .5, # Final
output_dir = "C:\output.las",
global_shift = F,
global_shift_type = "AUTO",
filter_sf = F,
filter_value = c(0, 1),
c_export_fmt = "LAS",
c_ext= "las", silent = T,
no_timestamp = T) {

if(class(file) != "character" ) {
stop(str_c("Wrong file is a ", class(file),", must be a character"))
}

if(class(filter_sf) != "logical" ) {
stop("Wrong filter_sf must be a logical value - TRUE or FALSE")
}

if(length(filter_value) != 2) {
stop("Wrong filter_value must be vector with 2 values for minimum and maximum treshold in filter_sf tool")
}

if(class(global_shift) != "logical" ) {
stop("Wrong global_shift must be a logical value - TRUE or FALSE")
}

if(!(global_shift_type %in% c("AUTO", "FIRST")) && (class(global_shift_type) != "numeric" || length(global_shift_type) != 3)) {
stop("Wrong global_shift_type must be character value : 'AUTO' or 'FIRST' (avaliable since CC v.2.11) or numeric vector with 3 values for dimensions x,y,z")
}

if(class(no_timestamp) != "logical" ) {
stop("Wrong no_timestamp must be a logical value - TRUE or FALSE")
}

if(!c_export_fmt %in% c("LAS", "LAZ")) {
stop(paste('Wrong c_export_fmt: ',c_export_fmt,' (',class(c_export_fmt),'), must be character value, one of these : LAS, LAZ', sep = ""))
}

if((!c_ext %in% c("LAS", "LAZ") && !c_ext %in% tolower(c("LAS","LAZ")))) {
stop(str_c('Wrong c_ext: ', c_ext,' (',class(c_ext),'),
must be character value, one of these :
LAS, LAZ
or
las, laz'))

}

auto_save_off = "-AUTO_SAVE OFF"

s1 = ifelse(silent, "-SILENT", "")

s2 = ifelse(no_timestamp, "-NO_TIMESTAMP", "")

s3 = str_c("-C_EXPORT_FMT", c_export_fmt, "-EXT", c_ext, sep = " ")

s4 = ifelse(global_shift, str_c("-O -GLOBAL_SHIFT", global_shift_type, file, sep = " "), str_c("-O", file, sep = " "))

s5 = str_c("-TREEISO", "-K1", K1, "-LAMBDA1", L1, "-DECIMATE_RESOLUTION1", DEC_R1,
"-K2", K2, "-LAMBDA2", L2, "-MAX_GAP", MAX_GAP, "-DECIMATE_RESOLUTION2", DEC_R2,
"-VERTICAL_OVERLAP_WEIGHT", VER_O_W, "-RHO", RHO, sep = " ")

s6 = ifelse(filter_sf, str_c("-FILTER_SF", filter_value[1], filter_value[2], sep = " "), "")

s7 = str_c("-SAVE_CLOUDS" , "FILE", output_dir, sep = " ")

cc_function = str_c(s1, s2, auto_save_off, s3, s4, s5, s6, s7, sep = " ")

return(cc_function)
}

#>>>>>>>
#> Example
CC_dir = "C:/CloudCompare/CloudCompare.exe" # path to CC.exe
indirdat = ".../A_01.las" # path to input las/laz file
outdir = ".../tryal" # path to output directory

#> Run
CC(cc_TREEiso(indirdat,
K1 = 10, L1 = 1, DEC_R1 = .1,
K2 = 20, L2 = 20, MAX_GAP = .2, DEC_R2 = .1,
VER_O_W = .5, RHO = .5
output_dir = str_c(outdir, "/tryal.las")))

#> Check
ltryal <- readLAS(str_c(outdir, "/tryal.las"))
ltryal@data$intermediate_segs %>% unique()
ltryal@data$final_segs %>% unique()
plot(ltryal, color = "final_segs", bg = "grey45", size = 3)
#>>>>>> END <<<<<<#