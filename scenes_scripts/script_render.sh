#!/bin/bash

#Check for dir, if not found create it 
output_name_file="Output/"
[ ! -d $output_name_file ] && mkdir $output_name_file

#Executable's location and name
exe_file_location="../build/Debug/"
exe_file_name="BunnyKiller.exe"

#Location of the scene to render
scene_files_location="../scenes/bunny/"

#Renderer parameters:
camera_samples=1
bdpt_samples=1
max_bounces=1

height=600
width=600

time_res=900
toffset=0
exp=0.025

camera_pos='1.6 0.5 -1.6'
camera_at='-0.5 0.0 1.5'

camera_fov=60

u_a='0.1'
u_s='0.1'

#Output file name
name='bunny_image'  
 
./$exe_file_location$exe_file_name \
-steady-state \
-log-name $output_name_file$name'_log.txt' \
-camera-spp $camera_samples -camera-fov $camera_fov \
-bidirectional-path-tracing $bdpt_samples \
-max-nb-bounces $max_bounces -steady-state \
-film-name $output_name_file$name -film-size-x $width -film-size-y $height \
-film-size-t $time_res -film-exposure $exp -film-offset $toffset \
-film-store-depth -film-store-normals -film-store-positions \
-camera-position $camera_pos -camera-focus $camera_at \
-point-light-source 0.0 1.0 -0.5 10.0 10.0 10.0 \
-name-mesh $scene_files_location'bunny.obj' -lambertian 0.5 0.0 0.5 \
-name-mesh $scene_files_location'wall.obj'  -lambertian 1.0 0.0 0.0 \
-name-mesh $scene_files_location'wall2.obj' -lambertian 0.7 0.7 0.7 \
-name-mesh $scene_files_location'floor.obj' -lambertian 0.7 0.7 0.7

