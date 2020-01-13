BunnyKiller Transient Renderer
==============================

BunnyKiller Transient Renderer is a small physically based renderer implementing  
the rendering algorithms developed for transient rendering described in:  

> "A Framework for Transient Rendering"  
> Adrian Jarabo, Julio Marco, Adolfo Muñoz, Raul Buisan, Wojciech Jarosz, Diego Gutierrez  
> ACM Transactions on Graphics, Vol.33(6), SIGGRAPH Asia 2014  

Overview
--------

The code compiles to a command line program that can render HDR images of scenes  
stored as obj files. The format of the command line arguments is explained in  
INPUT.md.

* The code uses the CImg library for output HDR images. Thus, while no  
  additional linking or dll is needed for output HDR format, other formats can  
  be used by using additional dlls.
* The code is designed with educational purposes, and therefore is not optimized  
  for speed, but to be the most comprehensive and extendible as possible.
* We will try to update the code with additional features and new goodies  
  periodically; unfortunately, this periodicity might not be as good as we'd  
  like. In any case, stay tuned to the
  [project page](http://giga.cps.unizar.es/~ajarabo/pubs/transientSIGA14/code).

License
-------

This piece of code is licensed under the MIT License (See LICENSE). Unless you  
really care, you can safely ignore all licenses in the source code.

Use it as you wish, but please reference the original publication:

    @article{Jarabo2014transient,
        author = {Jarabo, Adrian and Marco, Julio and
        	{Mu\~{n}oz}, Adolfo and Buisan, Raul and
        	Jarosz, Wojciech and Gutierrez, Diego},
        title = {A Framework for Transient Rendering},
        journal = {ACM Transactions on Graphics (SIGGRAPH Asia 2014)},
        volume = {33},
        number = {6},
        year = {2014},
    }

Contact
-------

 If you have any questions and/or suggestions for improving explanations or new  
 features, feel free to contact Adrian Jarabo <ajarabo@unizar.es>, or Julio  
 Marco <juliom@unizar.es>. Also, if you have coded any additional stuff that  
 you'd like us to include in the main project, we'll be happy to do it.

