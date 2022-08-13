# KalmanFilter :fire: :jack_o_lantern:

I found the tutorial of building this Kalman Filter on https://www.bzarg.com/p/how-a-kalman-filter-works-in-pictures/. I was unsure what several of the inputs would be so I followed my own way of calculating the covariance matrices and scaling matrices and yielded okay results.<br/>

The way this Kalman Filter works is it gathers data on points (x, y) and plugs them into a prediction equation. It is then compared to the actual results of the sensors (which I use random numbers for). There is a python file included and may be directly compiled from go.sh (please first use chmod u+x). Additionally you will need to download Matplotlib in order to plot in python, along with numpy and pandas. Below is the image of some of my results. For some reason the filter does well in predicting the X variable but not so much the Y variable. I would really appreciate some feedback on this error!

![alt](https://github.com/marscolony2040/KalmanFilter/blob/main/img/BB.png)<br/>

Pictured below is the math model used in this program which is found on the source url posted earlier
![alt](https://github.com/marscolony2040/KalmanFilter/blob/main/img/C.png)<br/>

## Running :cyclone:
In order to run this, you must have numpy, pandas, and matplotlib installed to access the 'kalman.py' file which plots your results. You do not need to install any additional Java libraries, because the only libraries imported are the default ones.<br/>

```sh
// Matrix algebra class
> javac Matrix.java

// Java statistical pack
> javac StatPack.java

// Kalman Filter
> javac kalman.java
> java kalman

// Plotting Kalman Results
> python kalman.py
```
