/* 
 * trim cube files to the minimal cube including all data points
 * beyond a given threshold. atom coordinates are preserved.
 *
 * optionally normalize a cube file for CPMD electrostatic potential 
 * calculations: calculate the average potential by summing over the 
 * grid points and subtract the result from all points.
 *
 * optionally select the phase (can be '+'/'p', '-'/'m'/'n' or '0'/'a[uto]')
 * and make the resulting data positive.
 *
 * Copyright (c) 2004-2005 by <axel.kohlmeyer@theochem.ruhr-uni-bochum.de>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int print_help(int retval, const char *opt)
{
    fprintf(stderr,"\ntrimcube: cut out a part of a cubefile.\n\n");

    if(retval==1) fprintf(stderr,"unknown option: '%s'\n\n", opt);
    if(retval==2) fprintf(stderr,"missing argument to option: '%s'\n\n", opt);
    if(retval==3) fprintf(stderr,"unknown argument to option: '%s'\n\n", opt);

    fprintf(stderr,"usage: trimcube [options] <input cube> [<output cube>]\n\n");
    fprintf(stderr,"available options:\n");
    fprintf(stderr,"-h[elp]          \tprint this info\n");
    fprintf(stderr,"-t[hresh] <value>\tset trim threshold to <value>, (default 0.005)\n");
    fprintf(stderr,"-n[orm]          \tnormalize so that the integral over the cube is zero.\n");
    fprintf(stderr,"-p[hase] <value> \tselect the phase. <value> can be '-','+', or '0'\n");
    fprintf(stderr,"\nuse '-' to read from standard input\n");
    fprintf(stderr,"program writes to standard output if not output file is given\n");
    fprintf(stderr,"\n");

    return retval;
}


int main(int argc, char **argv)
{
    int i, nra, nrb, nrc, nrdat, nat, *idx;
    int ia, ib, ic, mina, minb, minc, maxa, maxb, maxc;
    int lold, lnew, norm, phase;
    char titlea[255], titleb[255], *infile, *outfile;
    float origin[3], a[3], b[3], c[3], *data;
    float thresh, *q, *x, *y, *z, normdata, psum, nsum;
    FILE *fpin, *fpout;
    
    /* initialize variables */
    nat=nra=nrb=nrc=mina=minb=minc=norm=phase=0;
    thresh=0.005; psum=nsum=normdata=0.0;
    outfile=infile="-";
    fpin=fpout=NULL;
    
    /* parse arguments */
    i=1;
    while (i<argc) {
        if (0 == strncmp("-h", argv[i], 2)) return print_help(0, "");
        if (0 == strncmp("-n", argv[i], 2)) { norm=1; ++i; continue; }
        if (0 == strncmp("-t", argv[i], 2)) {
            ++i; 
            if (i>argc-1) return print_help(2, "-t");
            thresh=(float)atof(argv[i]);
            fprintf(stderr,"setting threshold to %f\n", thresh);
            ++i; continue;
        }
        if (0 == strncmp("-p", argv[i], 2)) {
            ++i;
            if (i>argc-1) return print_help(2, "-p");

            switch (argv[i][0]) {
              case 'a': /* fallthrough */
              case '0':
                  phase=999;
                  break;
                  
              case '-': /* fallthrough */
              case 'm':
              case 'n':
                  phase=-1;
                  break;

              case '+': /* fallthrough */
              case 'p':
                  phase=1;
                  break;

              default:
                  return print_help(3, "-p");
            }
            fprintf(stderr,"phase selection: %s\n", 
                    (phase<0) ? "negative" : (
                        (phase>1) ? "auto" : "positive") );

            ++i; continue;
        }
        

        /* stdin as file name */
        if (0 == strcmp("-", argv[i])) break;

        /* unknown option */
        if (0 == strncmp("-", argv[i], 1)) return print_help(1,argv[i]);

        /* filename for input cube */
        break;
    }

    /* we need at least an input file name ('-' for stdin). */
    if ((argc-i)<1) return print_help(0,"");
    infile=argv[i];
    if ((argc-i)>1) outfile=argv[i+1];

    /* open/assign input file pointer */
    if (0 == strcmp("-", infile)) {
        fpin = stdin;
    } else {
        fpin = fopen(infile, "r");
        if (!fpin) { 
            perror("could not open input cube"); 
            fprintf(stderr,"while trying to open '%s'\n", infile);
            return 10;
        }
    }
    
    /* read title strings */
    fgets(titlea, 255, fpin);
    fgets(titleb, 255, fpin);
    fscanf(fpin, "%d%f%f%f", &nat, origin, origin+1, origin+2);
    
    if (nat < 0) {
        fprintf(stderr,"trimcube: cannot handle orbital cube files\n");
        return 15;
    }
    fscanf(fpin, "%d%f%f%f", &nra, a, a+1, a+2);
    fscanf(fpin, "%d%f%f%f", &nrb, b, b+1, b+2);
    fscanf(fpin, "%d%f%f%f", &nrc, c, c+1, c+2);
    
    idx = (int *)malloc(nat*sizeof(int));
    q = (float *)malloc(nat*sizeof(float));
    x = (float *)malloc(nat*sizeof(float));
    y = (float *)malloc(nat*sizeof(float));
    z = (float *)malloc(nat*sizeof(float));

    for(i=0; i<nat; ++i) {
        fscanf(fpin, "%d%f%f%f%f", idx+i, q+i, x+i, y+i, z+i);
    }

    maxa=nra-1;
    maxb=nrb-1;
    maxc=nrc-1;
    nrdat = nra*nrb*nrc;
    data = (float *)malloc(nrdat*sizeof(float));
    for(i=0; i<nrdat; ++i) {
        if (1 != fscanf(fpin, "%f", data+i)) {
            fprintf(stderr,"trimcube: data read error on input cube\n");
            return 19;
        }
        normdata += data[i];

        /* for phase autodetect */
        if (data[i] > 0.0) {
            psum += data[i];
        } else {
            nsum += data[i];
        }
    }
    fclose(fpin);

    /* autodetect phase if not yet set.*/
    if (phase > 1) {
        if (psum+nsum > 0.0) {
            phase = 1;
        } else {
            phase = -1;
        }
        fprintf(stderr,"automatic phase detection gives: %s\n", 
                    (phase<0) ? "negative" : "positive" );
    }

    if (norm) {
        normdata = normdata / (float) nrdat;
        fprintf(stderr,"trimcube: subtracting average of %12.8f\n", normdata);
    } else {
        normdata = 0.0;
    }
    /* search for min/max region */
    if (thresh > 0.0) {
        maxa=maxb=maxc=-1;
        mina=nra;
        minb=nrb;
        minc=nrc;
        for (ia=0; ia<nra; ++ia) {
            for (ib=0; ib<nrb; ++ib) {
                for (ic=0; ic<nrc; ++ic) {

                    float d;

                    i=(ia*nrb + ib)*nrc + ic;
                    d = data[i];

                    /* manipulate data for phase selection */
                    if ( (phase>0) && (d<0.0) ){ d=0.0; }
                    if (phase<0) { d = (d>0.0) ? 0.0 : -d; }
                    if (phase!=0) data[i]=d;
                    
                    /* detect box with data above threshold */
                    if (((d > 0.0) && (d > thresh))
                        || ((d < 0.0) && (d < -thresh))) {
                        if (ia<mina) mina=ia;
                        if (ib<minb) minb=ib;
                        if (ic<minc) minc=ic;
                        if (ia>maxa) maxa=ia;
                        if (ib>maxb) maxb=ib;
                        if (ic>maxc) maxc=ic;
                    }
                }
            }
        }
    }
    
    /* sanity check */
    if (maxa < 0) {
	 fprintf(stderr, "trimming threshold is too high. " 
			 "refusing to write an empty cube file.\n");
	 return 99;
    }

    /* open/assign output file pointer */
    if (0 == strcmp("-", outfile)) {
        fpout = stdout;
    } else {
        fpout = fopen(outfile, "w");
        if (!fpout) { 
            perror("could not open output cube"); 
            fprintf(stderr,"while trying to open '%s'\n", outfile);
            return 20;
        }
    }

    /* write resulting cube. start with header. */
    fputs(titlea, fpout);
    fputs(titleb, fpout);

    origin[0] += mina*a[0] + minb*b[0] + minc*c[0];
    origin[1] += mina*a[1] + minb*b[1] + minc*c[1];
    origin[2] += mina*a[2] + minb*b[2] + minc*c[2];
    fprintf(fpout,"%5d%12.6f%12.6f%12.6f\n",nat,origin[0],origin[1],origin[2]);

    fprintf(fpout,"%5d%12.6f%12.6f%12.6f\n",maxa-mina+1,a[0],a[1],a[2]);
    fprintf(fpout,"%5d%12.6f%12.6f%12.6f\n",maxb-minb+1,b[0],b[1],b[2]);
    fprintf(fpout,"%5d%12.6f%12.6f%12.6f\n",maxc-minc+1,c[0],c[1],c[2]);
    for (i=0; i<nat; ++i) {
        fprintf(fpout,"%5d%12.6f%12.6f%12.6f%12.6f\n",idx[i],q[i],x[i],y[i],z[i]);
    }

    /* write remaining cube data */
    for (ia=mina; ia<=maxa; ++ia) {
        for (ib=minb; ib<=maxb; ++ib) {
            i=0;
            for (ic=minc; ic<=maxc; ++ic) {
                fprintf(fpout,"%13.5E", data[(ia*nrb + ib)*nrc + ic]
                        - normdata);
                ++i;
                if(i > 5) { fprintf(fpout,"\n"); i=0; }
            }
            if(i>0) fprintf(fpout,"\n");
        }
    }
    fclose(fpout);

    /* calculate compression (assuming unix CR/LF conventions). */
    lold = lnew = strlen(titlea)+strlen(titleb)+168+nat*54;
    lold += (nra*nrb*nrc)*13+((nrc/6+1)*nra*nrb);
    lnew += (maxa-mina+1)*(maxb-minb+1)*(maxc-minc+1)*13
        + (maxa-mina+1)*(maxb-minb+1)*((maxc-minc+1)/6+1);
    fprintf(stderr,"cube file reduced from %d to %d bytes. ratio: %3.1f:1\n",
           lold, lnew, ((double) lold)/((double) lnew));
    return 0;
}
