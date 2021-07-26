from oct2py import octave
octave.addpath(octave.genpath('/home/shiven/Desktop/GSoC-INCF/cdme/CDMentropy'))
octave.run('makeMex.m')
octave.run('try_me.m')