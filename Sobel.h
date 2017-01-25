//
//  Sobel.h
//  Proc
//
//  Created by Senda LI on 19/10/2016.
//  Copyright © 2016 IAOMW. All rights reserved.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/time.h> 
#include <sys/types.h>

#include "cv.h"
#include "opencv2/opencv.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;

Mat sobel_opt(Mat input);