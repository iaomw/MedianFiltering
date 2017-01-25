//
//  Proc.cpp
//  Proc
//
//  Created by Senda LI on 19/10/2016.
//  Copyright © 2016 IAOMW. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h> 

/* Libraries OPENCV */
#include "highgui.h"

#include "cv.h"
#include "opencv2/opencv.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include "Sobel.h"
#include "Median.h"

#ifndef max
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#endif
 
using namespace std;
using namespace cv;
using std::cout;

/*--------------- MAIN FUNCTION ---------------*/
int main () {
// préparation d'acquisition de flux vidéo
  VideoCapture cap(0);
  if(!cap.isOpened()){
    cout << "Errore"; return -1;
  }

// déclaration des variables 
// Mat structure contenant l'image
  Mat3b frame; Mat frame_gray;

    Mat grad_x, abs_grad_x;
    Mat grad_y, abs_grad_y;

  Mat frame_median, frame_sobel;
  Mat buffer_median, buffer_sobel;

  int ddepth = CV_16S;
  int scale = 1;
  int delta = 0;	
  char key = '0';

  cvNamedWindow("CV Median", WINDOW_AUTOSIZE);
  cvNamedWindow("CV Sobel", WINDOW_AUTOSIZE);
  cvNamedWindow("MY Median", WINDOW_AUTOSIZE);
  cvNamedWindow("MY Sobel", WINDOW_AUTOSIZE);

  cvMoveWindow("CV Median", 10, 30);
  cvMoveWindow("CV Sobel", 800, 30);
  cvMoveWindow("MY Median", 10, 500);
  cvMoveWindow("MY Sobel", 800, 500);
  
  while(key!='q'){
  // acquisition d'une trame video - librairie OpenCV
    cap.read(frame);

  //conversion en niveau de gris - librairie OpenCV
    cvtColor(frame, frame_gray, CV_BGR2GRAY);

    long long elapsed;
    struct timeval t_before, t_after;

    for (int i=3; i<32; i+=2) {

      medianBlur(frame_gray, frame_median, i);
      gettimeofday(&t_before, 0); 
      medianBlur(frame_gray, frame_median, i);
      gettimeofday(&t_after, 0); 

      elapsed = (t_after.tv_sec-t_before.tv_sec)*1000000LL + t_after.tv_usec-t_before.tv_usec;
      printf("Median CV elapsed for size %d is %lld\n", i, elapsed);

      buffer_median = median_s(frame_gray, i);
      gettimeofday(&t_before, 0); 
      buffer_median = median_s(frame_gray, i);
      gettimeofday(&t_after, 0); 

      elapsed = (t_after.tv_sec-t_before.tv_sec)*1000000LL + t_after.tv_usec-t_before.tv_usec;
      printf("MY Median S elapsed for size %d is %lld\n", i, elapsed);

      buffer_median = median_x(frame_gray, i);
      gettimeofday(&t_before, 0); 
      buffer_median = median_x(frame_gray, i);
      gettimeofday(&t_after, 0); 

      elapsed = (t_after.tv_sec-t_before.tv_sec)*1000000LL + t_after.tv_usec-t_before.tv_usec;
      printf("MY Median X elapsed for size %d is %lld\n", i, elapsed);

      frame_median = median_c(frame_gray, i);
      gettimeofday(&t_before, 0); 
      frame_median = median_c(frame_gray, i);
      gettimeofday(&t_after, 0); 

      elapsed = (t_after.tv_sec-t_before.tv_sec)*1000000LL + t_after.tv_usec-t_before.tv_usec;
      printf("MY Median C elapsed for size %d is %lld\n", i, elapsed);
    }

    gettimeofday(&t_before, 0); 
    medianBlur(frame_gray, frame_median, 11);
    gettimeofday(&t_after, 0); 
    elapsed = (t_after.tv_sec-t_before.tv_sec)*1000000LL + t_after.tv_usec-t_before.tv_usec;
    printf("elapsed for median CV %lld\n", elapsed);

    gettimeofday(&t_before, 0); 
    buffer_median = median_c(frame_gray, 11);
    gettimeofday(&t_after, 0); 
    elapsed = (t_after.tv_sec-t_before.tv_sec)*1000000LL + t_after.tv_usec-t_before.tv_usec;
    printf("elapsed for MY median C is %lld\n", elapsed);

    //calcul du gradient- librairie OpenCV
    gettimeofday(&t_before, 0);
    /// Gradient Y
    Sobel( frame_median, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
    convertScaleAbs( grad_x, abs_grad_x );
    /// Gradient Y
    Sobel( frame_median, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
    convertScaleAbs( grad_y, abs_grad_y );
    /// Total Gradient (approximate)
    addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, frame_sobel ); 
    gettimeofday(&t_after, 0); 
    elapsed = (t_after.tv_sec-t_before.tv_sec)*1000000LL + t_after.tv_usec-t_before.tv_usec;
    printf("elapsed for CV Sobel is %lld\n", elapsed);
    gettimeofday(&t_before, 0);

    gettimeofday(&t_before, 0);
    buffer_sobel = sobel_opt(frame_median);
    gettimeofday(&t_after, 0); 
    elapsed = (t_after.tv_sec-t_before.tv_sec)*1000000LL + t_after.tv_usec-t_before.tv_usec;
    printf("elapsed for MY Sobel is %lld\n", elapsed);
    
	// visualisation
	// taille d'image réduite pour meuilleure disposition sur écran
  // resize(frame, frame, Size(), 0.5, 0.5);
  // resize(frame_gray, frame_gray, Size(), 0.5, 0.5);
  // resize(grad, grad, Size(), 0.5, 0.5);

    imshow("CV Median",frame_median);
    imshow("CV Sobel", frame_sobel); 
    imshow("MY Median", buffer_median);  
    imshow("MY Sobel", buffer_sobel);

    //break;
    
    key=waitKey(5);
  }
}