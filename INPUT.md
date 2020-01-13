Bunnykiller input arguments
===========================

In general the system is unitless, everything is defined in terms of a  
'distance' unit, which is assumed to be meters by convention. The transient  
streak images capture light path length.

Arguments format
----------------

The renderer accepts arguments of the kind '-argument value value...'.

A value can be of the type:

* &lt;value: float\>: a floating point number separated by a '.' (e.g. 1.5).

* &lt;value: int\>: a integer number (e.g. 10).

* &lt;value: string\>: a text string, with properly string limiters (e.g. 'name'
or 'a long name').

*Warning*: The program expects all the input arguments to be correctly formatted  
and fails *silently* otherwise.

*Warning*: In case a parameter is specified twice or more, the last declaration  
will overwrite any previous one.

Camera arguments
----------------

* camera-position &lt;x: float\> &lt;y: float\> &lt;z: float\>

    x y z coords. of the camera position.

* camera-focus &lt;x: float\> &lt;y: float\> &lt;z: float\>

    x y z coords. of the point at which the camera is looking.

* camera-up &lt;x: float\> &lt;y: float\> &lt;z: float\>

    x y z coords. of the up direction of the camera plane.

* camera-spp &lt;spp: int\>

    Number of *squared* samples per pixel. The final sample count of camera  
    plane samples will be spp*spp;

* camera-fov &lt;angle: float\>

    Camera field of view (in degrees).

* camera-view-plane &lt;dist: float\>

    Distance from camera to view plane. Unless trying to imitate weird rendering  
    systems, better left unchanged.

Film arguments
--------------

By default the images are saved in Radiance format (.hdr)

* film-name &lt;name: string\>

    Name of the output file (without file extension). In steady state it will  
    result in a single 'name.hdr' file, in transient it will generate several  
    'name_XXXX.hdr' streak images, corresponding with scanlines along the height  
    of the image.

* film-size-x &lt;width: int\>

    Width of the image. If only this is indicated, the renderer assumes a square  
    image.

* film-size-y &lt;height: int\>

    Height of the image. 

* film-size-t &lt;length: int\>

    Length of the streak images.

* film-exposure &lt;exp: float\>

    Exposure per pixel. The distance the light can advance in each pixel of the  
    streak image. Every light path exceeding exp*length of the streak image will  
    be discarded.

* film-offset &lt;off: float\>

    Minimal path length that will be recorded by the streak images. It can be  
    interpreted like a time offset for the sensor.

* film-aspect-ratio &lt;ratio: float\>

    Forces a different aspect ratio for the image plane (by default image  
    width/height).

* film-scanline &lt;y: int\>

    Renders only the given scanline instead of the entire image

* film-single-pixel &lt;x: int\> &lt;y: int\>

    Renders only the given pixel instead of the entire image

* film-reveal

    Saves the image in plain text format (.csv) rather than in Radiance.

* -film-store-depth

    If set, a CImg RAW file containing the depth corresponding to each pixel will be generated.

* -film-store-positions

    If set, a CImg RAW file containing the position corresponding to each pixel will be generated.

* -film-store-normals

    If set, a CImg RAW file containing the geometric normal corresponding to each pixel will be  
    generated.

Transient arguments
-------------------

By default the renderer works in transient mode with no special time sampling.

* steady-state

    Forces steady state render. No streak images will be generated and there  
    will be no path length limit. The resulting intensity can be considerably  
    higher than the steady state image resulting from a transient render.

* transient-state

    Forces transient state render. 

* multibounce-streak &lt;N: int\>

    Stores up to N bounces in separate streak images, in the form  
    'name_bXX_XXXX.hdr'. Any radiance exceeding N bounces will be discarded.

* camera-unwarp

    Stores world-space distance rather than camera-space time.

Integrator arguments
--------------------

The renderer supports both bidirectional and progressive photon mapping  
integrators.

*Warning*: Selecting the integrator is a compile-time switch. Each integrator  
has its own parameters which are ignored otherwise.

* max-nb-bounces &lt;N: int\>

    Maximum number of bounces allowed in a single light path.

* scattering-level &lt;string\>

    Limits the number of participating media interaction allowed:

    * 'none': No media interactions.
    * 'single': A single media interaction.
    * 'multiple': Media interactions limited by the path length, ignoring the  
first one

    Not specifying this parameter will result in all participating media  
interactions being taken into account.

* time-sampling

    Forces time-aware path sampling heuristics (only useful in participating  
    media).

* standard-sampling

    Forces (default) path sampling techniques.

### Bidirectional Path tracing

* bidirectional-path-tracing &lt;N: int\>

    Generates N light paths for each camera sample.

### Progressive Photon Mapping

* photon-mapping &lt;N: int\> &lt;NN: int\> &lt;r: float\> &lt;M: int\>

    Traces N photons, uses the NN closest neighbours for density estimation for  
    each of the M eye paths and reduces the kernel radius by a factor of r in  
    each iteration.

Light arguments
---------------

* point-light-source &lt;x: float\> &lt;y: float\> &lt;z: float\> &lt;r: float\>
&lt;g: float\> &lt;b: float\>

    Omnidirectional point light source at position (x, y, z) that emits radiance  
    of intensity (r, g, b).

* point-light-source-polarized &lt;x: float\> &lt;y: float\> &lt;z: float\> 
&lt;I: float\> &lt;Q: float\> &lt;U: float\> &lt;V: float\>

    Point light source which can emit polarized monochromatic radiance  
    characterized by the Stokes vector [I Q U V]. Note that it emits in the   
    same reference frame as the sensor, so the camera must be specified  
    beforehand.

* cosine-light-source &lt;px: float\> &lt;py: float\> &lt;pz: float\> &lt;vx: float\> 
&lt;vy: float\> &lt;vz: float\> &lt;I: float\>

    Cosine-weighted light source at position (px, py, pz) that emits with  
    monochromatic radiance of intensity I in direction (vx, vy, vz).

* hemispherical-light-source &lt;px: float\> &lt;py: float\> &lt;pz: float\> 
&lt;vx: float\> &lt;vy: float\> &lt;vz: float\> &lt;I: float\>

    Light source at position (px, py, pz) that emits with monochromatic  
    radiance of intensity I in a hemisphere around direction (vx, vy, vz).

* spot-light-source &lt;px: float\> &lt;py: float\> &lt;pz: float\> &lt;vx: float\> 
&lt;vy: float\> &lt;vz: float\> &lt;ang: float\> &lt;I: float\>

    Light source at position (px, py, pz) that emits with monochromatic  
    radiance of intensity I in a cone of angle ang 8in degrees) around  
    direction (vx, vy, vz).

* rectangular-area-light-source &lt;px: float\> &lt;py: float\> &lt;pz: float\> 
&lt;vx: float\> &lt;vy: float\> &lt;vz: float\> &lt;width: float\> &lt;height: float\> 
&lt;I: float\>

    Rectangular area source, with lower corner (px, py, pz), a orientation  
    determined by its up direction (vx, vy, vz) and a extension width*height.

    By default it emits with monochromatic intensity I, but it can accept a  
    optional parameter '-emission-texture &lt;name: string>' so it emits with  
    a profile determined by the texture at file 'name'.

* directional-light-source &lt;x: float\> &lt;y: float\> &lt;z: float\> &lt;I: float\>

    A lights that emits uniformly over all the scene with intensity I coming  
    from direction (x, y, z).

* virtual-point-light-source

    A cosine-weighted light source. It can be specified by global coordinates  
    or relative to the camera:

    * &lt;cosine-light-source parameters\>

        It can be defined using the same parameter as a cosine light source.

    * -image-coordinates &lt;x: float\> &lt;y: float\> &lt;I: float\>

        It generates a light source  with monochrome intensity I at the first  
        intersection with the scene of a ray originating from the VPL emitter,  
        pointing towards the direction (x, y) over the image plane (defined  
        from -1 to 1, being 0 the center of the image). The light emites in the  
        direction of the normal at the hit point.

    This light source is mostly useful when defining DARPA lasers.

* virtual-point-light-emitter &lt;px: float\> &lt;py: float\> &lt;pz: float\>
&lt;dx: float\> &lt;dy: float\> &lt;dz: float\> &lt;ux: float\> &lt;uy: float\>
&lt;uz: float\>

    Emitter for virtual light sources, acts like a virtual camera at position  
    (px, py, pz) aiming at point (dx, dy, dz) and with up vector (ux, uy, uz).  
    It left unespecified, the actual camera will act as the virtual light  
    emitter.

    It can act like a perspective or a ortographic camera:

    * -ortographic &lt;w: float\> &lt;h: float\>

        Acts like a ortographic camera. The scale parameters w and h scale the  
        image plane coordinates of the virtual lights.

    * -perspective &lt;fov: float\>

        Acts like a perspective camera with field of view fov (in degrees).

Object arguments
----------------

Currently, the renderer uses Embree library internally for raytracing and only  
accepts triangle mesh geometry.

* name-mesh &lt;name: string\> &lt;material: material arguments\>

    Triangle mesh in Wavefront (.obj) format. How to define material parameters  
    is detailed in the next section. When no material parameter is provided, the  
    object uses the scene default material.

Material arguments
------------------

Certain materials can contain a participative media inside, whose parameters are  
defined in the next section

* lambertian &lt;r: float\> &lt;g: float\> &lt;b: float\>

    Lambertian material with reflectance (r, g, b).

* phong &lt;Kd: float\> &lt;Ks: float\> &lt;n: float\>

    Phong material with monochromatic diffuse reflectance Kd, monochromatic  
    specular reflectance Ks and specular parameter n.

* transparent &lt;r: float\> &lt;g: float\> &lt;b: float\>
&lt;medium: medium arguments\>

    Transparent interface which can contain a medium inside. When passing  
    through the interface light will experiment attenuation (r, g, b).

* dielectric &lt;r: float\> &lt;g: float\> &lt;b: float\> &lt;n1: float\> &lt;n0: float\>
&lt;medium: medium parameters\>

    Dielectric interface between a external medium with IOR n1 and a internal  
    medium with IOR n0.  If no medium is specified, light will experiment  
    attenuation (r, g, b).

* conductor &lt;n: float\> &lt;k: float\>

    Fresnel conductor material with real IOR n and imaginary IOR k.

* lambertian-textured &lt;rho: float\> &lt;tex: string\>

    Lambertian material material with absorption defined by the texture in file  
    'name', scaled by rho.

* two-sided &lt;material: material parameters\>

    Material adaptor which ignores geometric normals and acts as if the surface  
    is always facing the incident direction. Accepts any other material.

Medium arguments
----------------

Currently only homogeneous mediums are supported. They can use a variety of  
phase functions.

* homogeneous-medium &lt;ua_r: float\> &lt;ua_g: float\> &lt;ua_b: float\>
&lt;us_r: float\> &lt;us_g: float\> &lt;us_b: float\>
&lt;pf: phase function parameters\>

    Homogeneous medium with absorption (ua_r, ua_g, ua_b), scattering factor  
    (us_r, us_g, us_b) and a given phase function.

### Phase function arguments

All the phase functions have a optional '-pf-time &lt;float\>' parameter which  
adds a delay to every light-media interaction.

* pf-isotropic

    Radiance is remitted uniformly.

* pf-hg &lt;g: float\>

    Henyey-Greenstein phase function. The asymmetry parameter g interpolates  
    between backward and forward scattering.

* pf-rayleigh

    Rayleigh phase function

Scene-wide arguments
--------------------

* &lt;medium: medium-parameters\>

    When medium parameters without associated object are provided, all the  
    empty space in the scene is treated as filled with this participating media.

* lambertian-rho &lt;rho: float\>

    Scene default material, a lambertian material with monochromatic albedo  
    rho.

* log-name &lt;name: string\>

    Name of the file to save the render log.

* new-seed &lt;seed: int\>

    Initial seed for the internal RNG.

External files arguments
------------------------

* parse-args-file &lt;file: string\>

    Reads command-line arguments line by line from a external plain text file.  
    Mostly useful in Windows due to the command-line arguments limit.

* parse-xml-file &lt;file: string\>

    Loads a scene description from a Mitsuba-like XML format. The format is  
    documented in FORMAT.md
