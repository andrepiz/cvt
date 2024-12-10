# CVT
_Computer Vision Toolbox_

**Installation** 
Clone the repository by running the following git commands:

`git clone https://github.com/andrepiz/cvt`

Then, simply call _cvt_install()_ to add the toolbox to the PATH so you can call any CVT functions within your projects.

**Description** 
This toolbox contains many MATLAB functions to solve problems commonly found in the CV domain. Some examples are shown in the following pictures:

- Find the intersection of a Field of View (FOV) with conical or frustum shape with a sphere and compute the longitude/latitude limits at different camera poses. Note that when no intersections are found, the tangency limit circle is used to compute the limits.
![fov_intersection_scatter_pose](https://github.com/user-attachments/assets/bb503b7d-4a18-471f-a30d-f3efd7b7cdad)

- Given a sphere, camera and light locations, find the longitude/latitude masks of out-of-fov, non-observable and non-lit points. Non-observable points have a reflection angle larger than 90 degrees, non-lit points have a incidence angle larger than 90 degrees.
![image](https://github.com/user-attachments/assets/0dccc8f0-46ee-4ee0-9802-46c62ac49fe2)
![prova](https://github.com/user-attachments/assets/c5678d45-c68f-4d4f-9b42-2429cc69ae2e)

- Find the occlusions at limb or terminator regions of a sphere given its Digital Elevation Map (DEM). 
![occlusion_analysis](https://github.com/user-attachments/assets/41a32a74-57e5-47a8-8820-21e3d2c5d45d)

**Projects** 
This toolbox is currently used within the following projects:
- https://github.com/andrepiz/abram
