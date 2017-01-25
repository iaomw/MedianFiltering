//
//  Median.h
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

#include <stdio.h>
#include <pthread.h>

#include "cv.h"
#include "opencv2/opencv.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;
using namespace std;

Mat median_s(Mat input, int k_size);
Mat median_x(Mat input, int k_size);

Mat median_c(Mat input, int k_size);
