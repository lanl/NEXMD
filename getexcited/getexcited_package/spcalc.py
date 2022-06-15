#/usr/bin/python

'''

This function prepares input files for single-point calculations.

A general header called 'header' must be in the single-point
directory (e.g. singlepoint) and must have all inputs set except
for:

1) Initial nuclear coordinates and velocities (nucl_coord_veloc_flag)

In parentheses shown above, is the flag in 'header' that labels
this input.  This function finds this label and fills it in
accordingly.

NOTE: All NEXMD folders, inside the single-point directory, will
be deleted if this function is completely executed!

'''

import numpy as np
import random
import os
import sys
import shutil
import glob

def spcalc(header):

    print('Preparing input files for single-point calculations.')

    ## Directory names ##
    gsdir = input('Ground-state dynamics directory: ')
    if not os.path.exists(gsdir):
        print('Path %s does not exist.' % (gsdir))
        sys.exit()
    outdir = input('Output directory [e.g. singlepoint]: ')
    if not os.path.exists(outdir):
        print('Path %s does not exist.' % (outdir))
        sys.exit()

    ## Delete previous NEXMD folders ##
    NEXMDs = glob.glob('%s/NEXMD*/' % (outdir))
    NEXMDs.sort()
    if len(NEXMDs) != 0:
        contq = input('** WARNING ** All NEXMD folders inside %s will be deleted!\nContinue? Answer yes [1] or no [0]: ' % (outdir))
        if contq not in [1,0]:
            print('Answer must be 1 or 0.')
            sys.exit()
        if contq == 0:
            sys.exit()
    for NEXMD in NEXMDs:
        print('Deleting', '%s' % (NEXMD))
        shutil.rmtree(NEXMD)

    ## Information from header ##
    if not os.path.exists('%s/header' % (outdir)):
        print('Path %s/header does not exist.' % (outdir))
        sys.exit()
    header = header('%s/header' % (outdir))

    ## Check single-point calculation ##
    if header.n_class_steps != 0:
        print('User must change n_class_steps in %s/header to 0 for single-point calculations.' % (outdir))
        sys.exit()

    ## Find geometries ##
    if not os.path.exists('%s/coords.xyz' % (gsdir)):
        print('path %s/coords.xyz does not exist.' % (gsdir))
        sys.exit()
    datac = open('%s/coords.xyz' % (gsdir),'r')
    datac = datac.readlines()
    lenc = len(datac)
    if not os.path.exists('%s/velocity.out' % (gsdir)):
        print('Path %s/velocity.out does not exist.' % (gsdir))
        sys.exit()
    datav = open('%s/velocity.out' % (gsdir),'r')
    datav = datav.readlines()
    lenv = len(datav)
    ncoords = 0
    index = 0
    arrayc = np.array([])
    for line in datac:
        if 'time' in line:
            if ncoords == 0:
                tinit = np.float(line.split()[-1])
            else:
                time = np.float(line.split()[-1])
            ncoords += 1
            arrayc = np.append(arrayc,index)
        index += 1
    arrayc = np.append(arrayc,lenc + 1)
    arrayc = np.int_(arrayc)
    if ncoords == 0:
        print('No coordinates were found.')
        sys.exit()
    if ncoords == 1:
        print('Only initial coordinates, at %.2f fs, were found.' % (tinit))
        sys.exit()
    if ncoords > 1:
        index = 0
        arrayv = np.array([])
        for line in datav:
            if 'time' in line:
                arrayv = np.append(arrayv,index)
            index += 1
        arrayv = np.append(arrayv,lenv)
        arrayv = np.int_(arrayv)
    tinc = time/(ncoords - 1)
    print('A total of %d coordinates files, ranging from %.2f to %.2f fs in increments of %.2f fs, were found.' % (ncoords,tinit,time,tinc))

    ## Choose geometries ##
    #coords = input('Enter requested range of the ground-state sampling by coordinate files and the number of single-point calculations.\nInput an array of the form [start, end, number]: ')
    coords = list(map(int,input("\n'Enter requested range of the ground-state sampling by coordinate files and the number of single-point calculations.\nInput an array of the form [start, end, number]: ").strip().split()))
    if not isinstance(coords,list):
        print('Input must be an array of the form [start, end, number].\nFor example, [1, 1000, 500] requests 500 coordinate files sampled from 1 to 1000.')
        sys.exit()
    if len(coords) != 3:
        print('Input must be an array with three elements of the form [start, end, number].\nFor example, [1, 1000, 500] requests 500 coordinate files sampled from 1 to 1000.')
        sys.exit()
    index = 0
    for i in coords:
        if isinstance(i, int) == False:
            print('Element number %d of input array must be integer.\nUser inputted [%s, %s, %s], which is not allowed.' % (index + 1,coords[0],coords[1],coords[2]))
            sys.exit()
        if index in [0,1]:
            if i not in np.arange(ncoords):
                print('Element number %d of input array must be integer between 0 and %d.\nUser inputted [%d, %d, %d], which is not allowed.' % (index + 1,ncoords - 1,coords[0],coords[1],coords[2]))
                sys.exit()
        else:
            if i not in np.arange(1,ncoords + 1) and i not in -np.arange(1,ncoords + 1):
                print('Element number %d of input array must be integer between 1 and %d.\nUser inputted [%d, %d, %d], which is not allowed.' % (index + 1,ncoords,coords[0],coords[1],coords[2]))
                sys.exit()
        index += 1
    if coords[0] > coords[1]:
        print('second element of input array must be greater than first element.\nUser inputted [%d, %d, %d], which is not allowed.' % (coords[0],coords[1],coords[2]))
        sys.exit()
    if (coords[1] - coords[0]) + 1 < np.abs(coords[2]):
        print('Number of coordinate files requested must be less than or equal to the number of coordinate files in the sample.\nUser inputted [%d, %d, %d], which is not allowed.' % (coords[0],coords[1],coords[2]))
        sys.exit()
    if coords[2] < 0:
        ntraj = np.abs(coords[2])
        coordsq = input('You have requested %d randomly-selected coordinate files in the range %d to %d.\nContinue? Answer yes [1] or no [0]: ' % (ntraj,coords[0],coords[1]))
    else:
        interval = np.int(np.ceil((coords[1] - coords[0] + 1)/np.float(coords[2])))
        if interval*coords[2] > (coords[1] - coords[0] + 1):
            interval = interval - 1
            ntraj = len(np.arange(coords[0],coords[1] + 1,interval))
        else:
            ntraj = coords[2]
        coordsq = int(input('You have requested %d evenly-spaced coordinate files in the range %d to %d.\nContinue? Answer yes [1] or no [0]: ' % (ntraj,coords[0],coords[1])))
    if coordsq not in [1,0]:
        print('Answer must be 1 or 0.')
        sys.exit()
    if coordsq == 0:
        sys.exit()

    ## Split geometries ##
    split = int(input('Number of single-point calculations per NEXMD folder [e.g. 100]: '))
    if isinstance(split, int) == False:
        print('Number of single-point calculations per NEXMD folder must be integer.')
        sys.exit()
    if split < 0:
        print('Number of single-point calculations per NEXMD folder must be integer greater than zero.')
        sys.exit()
    dirsplit = split*np.arange(1,np.ceil(np.float(ntraj)/split) + 1)
    dirsplit[-1] = coords[1] + 1
    if coords[2] < 0:
        dirsplit = np.split(np.sort(random.sample(np.arange(coords[0],coords[1]+1),ntraj)),dirsplit)
    else:
        dirsplit = np.split(np.arange(coords[0],coords[1]+1,interval),dirsplit)

    ## Extract atomic numbers ##
    if not os.path.exists('%s/restart.out' % (gsdir)):
        print('Path %s/restart.out does not exist.' % (gsdir))
        sys.exit()
    anum = open('%s/restart.out' % (gsdir),'r')
    anum = anum.readlines()
    top = None
    bottom = None
    index = 0
    for line in anum:
        if '$COORD' in line:
            top = index
        if '$ENDCOORD' in line:
            bottom = index
            break
        index += 1
    if top == None or bottom == None:
        print('There is a problem with %s/restart.out.' % (gsdir))
        sys.exit()
    if isinstance(top, int) == True and isinstance(bottom, int) == True:
        anum = [ line.split()[0] for line in anum[top+1:bottom:1] ]
    else:
        print('There is a problem with %s/restart.out.' % (gsdir))
        sys.exit()

    ## Prepare single-point input files ##
    index = 0
    for NEXMD in np.arange(1,np.ceil(np.float(ntraj)/split)+1):
        os.makedirs('%s/NEXMD%d' % (outdir,NEXMD))
        dirlist = open('%s/NEXMD%d/dirlist' % (outdir,NEXMD),'w')
        for dir in dirsplit[index]:
            coords = datac[arrayc[dir]+1:arrayc[dir+1]-1:1]
            velocs = datav[arrayv[dir]+2:arrayv[dir+1]-1:1]
            os.makedirs('%s/NEXMD%d/%04d' % (outdir,NEXMD,dir))
            inputfile = open('%s/NEXMD%d/%04d/input.ceon' % (outdir,NEXMD,dir),'w')
            for line in header.file:
                if 'nucl_coord_veloc_flag' in line:
                    inputfile.write('&coord\n')
                    aindex = 0
                    for line in coords:
                        val = line.split()
                        inputfile.write('{:>6}  {:>12}  {:>12}  {:>12}'.format(anum[aindex],val[1],val[2],val[3]))
                        inputfile.write('\n')
                        aindex += 1
                    inputfile.write('&endcoord\n\n&veloc\n')
                    for line in velocs:
                        inputfile.write(line)
                    inputfile.write('&endveloc\n')
                else:
                    inputfile.write(line)
            dirlist.write('%04d' % (dir))
            print('%s/NEXMD%d/%04d' % (outdir,NEXMD,dir))
        dirlist.close()
        shutil.copyfile('%s/NEXMD%d/dirlist' % (outdir,NEXMD), '%s/NEXMD%d/dirlist1' % (outdir,NEXMD))
        index += 1
