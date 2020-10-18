## Digital Image Processing
# Project Proposal
## Project ID:: 32
## Project Title:: Single Image Motion Deblurring Using Transparency
Paper Link:: <http://jiaya.me/all_final_papers/motion_deblur_cvpr07.pdf> 

## GtiHub repo link
<https://github.com/Digital-Image-Processing-IIITH/project-proximacentauri> 

## Team Members
* Abhinaba Bala **(2020701001)**, abhinaba.bala@research.iiit.ac.in
* Neel Ashok Mishra **(2020701009)**, neel.mishra@research.iit.ac.in
* Ritam Basu **(2020701005)**, ritam.basu@research.iiit.ac.in
* Jigyasu Khandelwal **(2020702013)**, jigyasu.khandelwal@research.iiit.ac.in

## Main Goals of the project
1. Image blur is caused due to object motion or camera vibration when the shutter is pressed. 
   The relationship between the blur filter and the object boundary transparency has to be understood. 
2. An algorithm is to be developed to estimate the blur filter from the transparency values given a blurred image. 
   This formulation is to be used to estimate a filter for general transparency present in an image. 
3. This will alow us to deblur an image irrespective of the cause of the blur, which is a first of its kind method (when) formulated by the authors. 

## Problem Definition 
1. Motion blur is usually formulated as a linear image degradation process, given by 
    _I = L * f + n_, where _I_, _L_ and _n_ represent the degraded captured image, unblurred (or latent) image, and the additive noise respectively. _*_ is the convolution operator and _f_ is the unkown linear shift-invariant point spread function. 
2. Conventional _blind deconvolution_ methods estimate the blur filter from image intensities or gradients and deconvolve the blurred image. 
3. A shortcoming of this method is that it can't completely solve the probelm because the background may not undergo the same motion as the object, as the filter is defined only on the moving object. 
4. The relationship between image boundary transparency and the deblur filter is investigated and it is shown that a filter can be estimated from only the transparency values. 
5. An algorithm is developed to estimate the filter using a Maximum A Prosteriori (MAP) formulation with a suitable prior and likelihood on transparency. 
6. A formulation of general transparency is investigated and it is shown that the previous formulation is robust for any kind of motino blur. 
7. The algorithm is also very efficient since only image patches are used as input. 

## Algorithms
1. Iterative optimization method for our MAP approach to recover the motion blur filter using transparency.  
Employ conjugate gradient optimization and then Belief Propagation to estimate the filter
2. Deconvolution using Lucy-Richardson method.
3. Use user-drawn strokes to collect the foreground and background samples for object motion blur.

## Expected Results of the project
1. We look to obtain a general formulation to solve the problem of image blurring using transparency for _i._ object blur and/or _ii._ camera blur.
2. The accurate estimation of transparency values is very crucial to the success of our method. 

## Project Milestones and Expected Timeline
_Project Proposal Submission_: 18th October  
25th October: Analysis of Paper  
31st October: Set up the starter code and dependencies  
_Mid Evaluation_: 31st October  
7th November: Implement MAP approach to recover the motion blur filter using transparency and solve 2-D object motion blur.   
14th November: Implement generalized transparency for when the entire image is degraded.  
19th November: Complete the codebase and prepare presentation.  
_Final Evaluation_: 19th-25th November  

## Is there a dataset that we require?
No
