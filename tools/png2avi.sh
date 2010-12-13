#!/bin/bash
# from http://www.mplayerhq.hu/DOCS/HTML/en/menc-feat-enc-images.html

mencoder mf://*.png -mf w=800:h=600:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
