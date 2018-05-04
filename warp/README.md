# Installation
The installation instructions for Warp are [here](http://warp.lbl.gov/home/how-to-s/installation).

Make sure to follow these with your sandbox turned on (`source ~/.virtual_envs/warp/bin/activate`).  As you already should have numpy and Forthon, you need not follow those portions of the instructions.  Also, you need not do the complicated parallel installation instructions, and stop after installing pygist.  

# Testing
Next, go to the bottom of the installation instruction page and follow the instructions for running the tests.  When I ran the tests, I got one error in pyelem_test.py.  However, my warp runs fine.

Once you have warp installed, you can see how Warp runs in the following manner:  Let $INSTALLATION_DIR represent the path to the source directory from the installation instructions.  Create a new directory wherever you want to store files related to this course.  In that directory and with the warp sandbox turned on, run `python $INSTALLATION_DIR/warp/scripts/examples/Pierce_diode.py`.  You should see some output, and after a couple seconds, the shell prompt.  Look at the files in this directory.  You should see one called Pierce_diode.000.cgm.  You can look at this file using pygist, `gist Pierce_diode.000.cgm`.  A document will open in an X-window, and you can navigate forward and backward through the document with f and b.  Move forward, and you should eventually see some snapshots of a current (black dots) coming out of a Pierce Geometry diode with blue lines representing electric fields as the current is becoming steady state.  This is one of the plots you should see:

![Pierce diode stead state current](https://raw.githubusercontent.com/billyziege/USPAS_Summer_2018/master/warp/pierce_diode.png)
  
After you close the window, you may need to “Ctrl C” to get out of the gist prompt and return to the shell prompt.
