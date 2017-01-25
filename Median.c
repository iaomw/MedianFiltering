//
//  Median.c
//  Proc
//
//  Created by Senda LI on 19/10/2016.
//  Copyright © 2016 IAOMW. All rights reserved.
//

#include "Median.h"

void swap_u8(uint8_t* v, int i, int j) {

  uint8_t c;
  c = v[i];
  v[i] = v[j];
  v[j] = c;
}

uint8_t q_s(uint8_t* v, int left, int right, int p_m) {

  while(left<right) {
    // the last index which 
    // we are sure it's no bigger than selected value.
    int last = left;
    int mid = (left+right)/2;

    swap_u8(v, left, mid);

    for (int i=left+1; i<=right; i++) {
      if (v[i] < v[left]) { //Swap that just next to the last smaller value.
        swap_u8(v, last++, i);
      }
    }

    swap_u8(v, left, last);
    if (last == p_m) { return v[last]; }
    else if (last > p_m) { right = last-1; }
    else if (last < p_m) { left = last+1; }
  }

  return v[p_m];
} 

int* row_dict;
uint8_t* copys;
uint8_t* sample;

int* histogram;

static int kernal_length;

Mat median_x(Mat input, int k_size) {

  if (k_size<3||k_size%2==0) {return input;}

  int radius = (k_size-1)/2;
  int capacity = pow(k_size,2);
  int p_median = (capacity-1)/2;

  Mat output = input.clone();
  int c = input.cols;
  int r = input.rows;

  int should_change = (k_size!=kernal_length);

  if (sample!=NULL && should_change) {
    free(sample);sample = NULL;
    free(copys);copys = NULL;
    free(row_dict);
    row_dict=NULL;
  }

  if (histogram==NULL) {
    histogram = (int*)calloc(256, sizeof(int));
  }

  if (sample==NULL) {
    sample = (uint8_t*)calloc(capacity, sizeof(uint8_t));
    copys = (uint8_t*)calloc(capacity, sizeof(uint8_t));
    row_dict = (int*)calloc(k_size, sizeof(int));
  }

  int x_s = radius, x_e = c-1-radius, x_step = 1;
  int y_s = radius, y_e = r-1-radius, y_step = 1;

  for (int e=0; e<k_size; e++) {
      int shift = x_s-radius+(y_s-radius+e)*c;
      memcpy(&sample[1+e*k_size], &*(input.data+shift), (k_size-1)*sizeof(uint8_t));
  } // [1+e*k_size] reserve one position for first iteration.

  int remain = k_size;
  for (int x=x_s; x<=x_e; x+=x_step) {
    // There will be a problem if image rows count is smaller than 2 time the kernal height.
    // Because while switch, different row has different data position.
    for (int e=0; e<k_size; e++) {
      sample[e*k_size] = *(input.data+x+radius+c*row_dict[e]);
    }

    for (int y=y_s; (y_step>0?y<=y_e:y>=y_e); y+=y_step, remain++) {

      remain = remain % k_size;
      remain = (remain == 0 ?k_size:remain);
      int y_shift = y+y_step*radius;
      int shift = y_shift*c+x-radius;
      int row = (remain-1);
      row_dict[row] = y_shift;

      memcpy(&sample[row*k_size], &*(input.data+shift), k_size*sizeof(uint8_t));
      //memcpy(copys, sample, capacity*sizeof(uint8_t));
      //int me = q_s(copys, 0, capacity-1, p_median);
      for (int e=0; e<capacity; e++) {
        histogram[sample[e]]++;
      }

      // do Median
      int me, count = 0;
      for (me=0; (me<256 && count<p_median); me++) {
        if (histogram[me]==0) { continue; }
        count += histogram[me];
      }

      memset(histogram, 0, 256*sizeof(int));

      *(output.data+x+y*c) = me;
    }// y

    y_step=-y_step;
    int t = y_s;
    y_s = y_e;
    y_e = t;
  }// x
  
  return output;
}

Mat median_s(Mat input, int k_size) {

  if (k_size<3||k_size%2==0) {return input;}

  int radius = (k_size-1)/2;
  int capacity = pow(k_size,2);
  int p_median = (capacity-1)/2;

  Mat output = input.clone();
  int c = input.cols;
  int r = input.rows;

  int should_change = (k_size!=kernal_length);

  if (sample!=NULL && should_change) {
    free(sample);sample = NULL;
    free(copys);copys = NULL;
    free(row_dict);
    row_dict=NULL;
  }

  if (sample==NULL) {
    sample = (uint8_t*)calloc(capacity, sizeof(uint8_t));
  }

  int x_s = radius, x_e = c-1-radius, x_step = 1;
  int y_s = radius, y_e = r-1-radius, y_step = 1;

  for (int x=x_s; x<=x_e; x+=x_step) {
    for (int y=y_s; y<=y_e; y+=y_step) {
      int shift = x-radius+(y-radius)*c;
      for (int e=0; e<k_size; e++) {
         memcpy(&sample[e*k_size], &*(input.data+shift+e*c), k_size*sizeof(uint8_t));
      }

      int me = q_s(sample, 0, capacity-1, p_median);
      *(output.data+x+y*c) = me;
    }// y
  }// x
  
  return output;
}

static Mat cvin;
static Mat cvout;

static int last_k_size;
static int last_i_cols;
static int last_i_rows;

uint16_t** h_cols;

struct thread_args {
    int k_size;
    uint16_t* h_cent;
    int xs, xe, ys, ye;
    int radius, capacity, p_median;
};

void* median_c_part(void* args) {

  struct thread_args *params = (thread_args*)args;

  uint16_t* histogram_cent = params->h_cent;
  uint16_t** histogram_cols = h_cols;

  int k_size = params->k_size;
  int radius = params->radius;
  int capacity = params->capacity;
  int p_median = params->p_median;

  int xs = params->xs; int xe = params->xe;
  int ys = params->ys; int ye = params->ye;

  int c = cvin.cols; //int r = cvin.rows;

  // Should only for border;
  int x_s = xs+radius+1, x_e = xe-1-radius;
  int y_s = radius+1, y_e = ye-1-radius;

  //free(args);

  for(int x=xs; x<=xe; x++) {
    memset(histogram_cols[x], 0, 256*sizeof(uint16_t));
    for (int y=(y_s-radius-1); y<(y_s+radius); y++) {
      Scalar ad = cvin.at<uchar>(Point(x, y));
      histogram_cols[x][(int)ad.val[0]]++;
    }
  } // reset

  int x_step = 1, y_step = 1;
  for (int x=(x_s-radius-1); x<(x_s+radius); x++) {
          for (int e=0; e<256; e++) {
              histogram_cent[e]+=histogram_cols[x][e]; e++;
              histogram_cent[e]+=histogram_cols[x][e]; e++;
              histogram_cent[e]+=histogram_cols[x][e]; e++;
              histogram_cent[e]+=histogram_cols[x][e]; e++;
              histogram_cent[e]+=histogram_cols[x][e]; e++;
              histogram_cent[e]+=histogram_cols[x][e]; e++;
              histogram_cent[e]+=histogram_cols[x][e]; e++;
              histogram_cent[e]+=histogram_cols[x][e];
          }
  }

  for (int y=y_s; y<=y_e; y+=y_step) {
    int xs, xe, d_intake, d_waive;
    if (x_step > 0) {
      xs = x_s-radius-1;
      d_intake = radius;
      d_waive = radius+1;
    } else {
      xs = x_s-radius+1;
      d_intake = -radius;
      d_waive = -radius-1;
    } xe = xs+k_size;

    for (int e=xs; e<xe; e++) {

      Scalar de = cvin.at<uchar>(Point(e, y-radius-1));
      histogram_cols[e][(int)de.val[0]]--;
      histogram_cent[(int)de.val[0]]--;
      Scalar ad = cvin.at<uchar>(Point(e, y+radius));
      histogram_cols[e][(int)ad.val[0]]++;
      histogram_cent[(int)ad.val[0]]++;
    }

    for (int x=x_s; (x_step>0?x<=(x_e+1):x>=(x_e-1)); x+=x_step) {
        int intake = x+d_intake, waive = x-d_waive;

        Scalar de = cvin.at<uchar>(Point(intake, y-radius-1));
        histogram_cols[intake][(int)de.val[0]]--;
        Scalar ad = cvin.at<uchar>(Point(intake, y+radius));
        histogram_cols[intake][(int)ad.val[0]]++;

              for (int e=0; e<256; e++) {
                histogram_cent[e] -= histogram_cols[waive][e];
                histogram_cent[e] += histogram_cols[intake][e]; e++;
                histogram_cent[e] -= histogram_cols[waive][e];
                histogram_cent[e] += histogram_cols[intake][e]; e++;
                histogram_cent[e] -= histogram_cols[waive][e];
                histogram_cent[e] += histogram_cols[intake][e]; e++;
                histogram_cent[e] -= histogram_cols[waive][e];
                histogram_cent[e] += histogram_cols[intake][e];
              } // shift histogram_cent

      // do Median
      int me, count = 0;
      for (me=0; (me<256 && count<p_median); me++) {
        if (histogram_cent[me]==0) { continue; }
        count += histogram_cent[me];
      }

      //cvout.at<uchar>(y*c+x) = me;
      cvout.at<uchar>(Point(x, y)) = me;
    }// x

    x_step=-x_step;
    int t = x_s;
    x_s = x_e;
    x_e = t;
  }// y 

  return NULL;  
}

static struct thread_args** t_info;
Mat median_c(Mat input, int k_size) {

  if (k_size<3||k_size%2==0) {return input;}

  int radius = (k_size-1)/2;
  int capacity = pow(k_size, 2);
  int p_median = (capacity+1)/2;

  cvout = input.clone(); cvin = input;
  int c = input.cols, r = input.rows;

  // for (int i=0; i<c; i++) {
  //   for (int j=0; j<r; j++) {
  //     cvout.at<uchar>(Point(i, j)) = 0;
  //   }
  // }

  int t_number = 4;

  if (t_info==NULL) {
    t_info = (struct thread_args**)malloc(t_number*sizeof(struct thread_args*));
    for (int i=0; i<t_number; ++i) {
      t_info[i] = (struct thread_args*)malloc(sizeof(struct thread_args));
      t_info[i]->h_cent = (uint16_t*)malloc(256*sizeof(uint16_t)); 
    }
  } 

  int should_change = c!=last_i_cols;

  last_k_size = k_size; 
  last_i_cols = c; 
  last_i_rows = r;

  if (should_change && h_cols!=NULL) {
    for(int i=0; i<last_i_cols; i++) {
      free(h_cols[i]);
      h_cols[i]=NULL;
    }
    free(h_cols);
    h_cols=NULL;
  } 

  if (h_cols==NULL) {
    h_cols = (uint16_t**)malloc(c*sizeof(uint16_t*));
    for(int x=0; x<c; x++) {
      h_cols[x] = (uint16_t*)malloc(256*sizeof(uint16_t));
    } // x
  } // NULL

  pthread_t* t_pool = (pthread_t*)malloc(t_number*sizeof(pthread_t));

  //int x_s = radius+1, x_e = c-2-radius;
  //int y_s = radius+1, y_e = r-2-radius;

  int x_s = 0, x_e = c-1;
  int y_s = 0, y_e = r-1;

  int piece = c/t_number-1;

  for (int i=0; i<t_number; i++, x_s+=1) {

    t_info[i]->k_size = k_size; 
    t_info[i]->radius = radius; 
    t_info[i]->capacity = capacity; 
    t_info[i]->p_median = p_median;

    memset(t_info[i]->h_cent, 0, 256*sizeof(uint16_t));

    t_info[i]->xs = x_s; 
        x_s+=piece;
    t_info[i]->xe = x_s;
    t_info[i]->ys = y_s; 
    t_info[i]->ye = y_e;

    pthread_create(&t_pool[i], NULL, median_c_part, (void*)t_info[i]);
  }

  for (int i=0; i<t_number; i++) { 
    pthread_join(t_pool[i], NULL);
  }

  return cvout;
}

//__m128i v_ce, v_co;
          // //for (int e=0; e<256; e+=16) {
          //     // v_ce = _mm_loadu_si128((__m128i*)(&histogram_cent[e]));
          //     // v_co = _mm_loadu_si128((__m128i*)(&histogram_cols[x][e]));
          //     // _mm_storeu_si128((__m128i*)(&histogram_cent[e]), _mm_adds_epu8(v_ce, v_co));
          //     // v_ce = _mm_loadu_si128((__m128i*)(&histogram_cent[e+16]));
          //     // v_co = _mm_loadu_si128((__m128i*)(&histogram_cols[x][e+16]));
          //     // _mm_storeu_si128((__m128i*)(&histogram_cent[e+16]), _mm_adds_epu8(v_ce, v_co));

 //__m128i v_ce, v_co, v_cc;
              // for (int e=0; e<256; e+=16) {
              //   v_ce = _mm_loadu_si128((__m128i*)(&histogram_cent[e]));
              //   v_co = _mm_loadu_si128((__m128i*)(&histogram_cols[waive][e]));
              //   v_cc = _mm_subs_epu8(v_ce, v_co);
              //   v_co = _mm_loadu_si128((__m128i*)(&histogram_cols[intake][e]));
              //   v_cc = _mm_adds_epu8(v_co, v_co);
              //   _mm_storeu_si128((__m128i*)(&histogram_cent[e]), v_cc);