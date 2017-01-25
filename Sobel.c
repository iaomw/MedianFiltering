//
//  Sobel.c
//  Proc
//
//  Created by Senda LI on 19/10/2016.
//  Copyright © 2016 IAOMW. All rights reserved.
//

#include "Sobel.h"

Mat sobel_opt(Mat input) {

  Mat output = input.clone(); 
  uchar *im_in = input.data;
  int n = input.cols;
  int m = input.rows;
  int s = n*m;

  int N,NW,NE, NE2,W2, W,E,E2,SW,SE,SE2,S;
  double gx, gx2;
  double gy, gy2;
  double root, root2;

  int i;
  for (i=1;i<(s-2);i+=2)
    { 
      NW = *(im_in+max(i-n-1,0));
      N = *(im_in+max(i-n,0));

      NE = *(im_in+max(i-n+1,0));     
      NE2 = *(im_in+max(i-n+2,0));
     
      W = *(im_in+i-1);
      W2 = *(im_in+i);
      E = *(im_in+i+1);
      E2 = *(im_in+i+2);
    
      SW = *(im_in+(i+n-1));
      S = *(im_in+(i+n));
      SE = *(im_in+(i+n+1));
      SE2 = *(im_in+(i+n+2)); 

      gx = NW + NE + 2*N - SW - SE - 2*S;
      gy = NE + SE + 2*E - NW - SW - 2*W;
      gx = abs(gx);
      gy = abs(gy);
      root = sqrt((gy*gy) + (gx*gx));
      *(output.data+i) = (uchar)root;
      //output.at<uchar>(i) = (uchar)root;
       
      gx2 = N + NE2 + 2*NE - S - SE2 - 2*SE; 
      gy2 = NE2 + SE2 + 2*E2 - N - S - 2*W2;
      gx2 = abs(gx2);
      gy2 = abs(gy2);
      root2 = sqrt((gy2*gy2) + (gx2*gx2));
      *(output.data+i+1) = (uchar)root2;
    }

    return output;
}