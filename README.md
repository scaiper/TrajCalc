TrajCalc
========

Trajectory Calculator for interplanetary flights with multiple gravity assists

This tool is designed for calculating flight plans. Given a sequence of planets, departure date range and maximum time of flight, it computes optimal flight plan in terms of dV or dV and time of flight simultaneously. In later case, result is a set of solutions represented as 2D plot, and user can pick one with reasonable trade off between dV and time of flight.

It has no GUI yet, you put flight parameters in config file and launch the application. Read comments in config.txt for details.

Results are output in two formats. Text in *.txt files. And plots in gnuplot format in *.plt files, you need gnuplot (http://www.gnuplot.info/) to visualize them.

Source code is available at GitHub under GPL. https://github.com/scaiper/TrajCalc 
