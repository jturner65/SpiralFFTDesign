package SpiralFFTDesPKG;


import java.util.*;
import org.apache.commons.math3.complex.*;
import org.apache.commons.math3.transform.*;

import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PImage;
//process and visualize fft of screen to use as tool 
public class myFFTCalcAndVis implements Runnable {
	public static FFTSpiralDsn pa;
	
	public Complex[] fftRes;
	public FastFourierTransformer fft;
	public int[] pxlArray;
	
	public int pxlAraSize, dispImgDim;
	public ArrayList<Integer> paddedZero;
	
	public PImage dispImage;
	
	public static int maxDispVals, minDispVals, max2ByteVals, min2ByteVals;
	
	public boolean finishedCalc;
	
	
	public myFFTCalcAndVis(FFTSpiralDsn _pa, PImage src) {
		pa = _pa;
		fft = new FastFourierTransformer(DftNormalization.STANDARD);
		int origPxlSize = src.pixels.length;
		int pwr = ((int) (Math.log(origPxlSize)/Math.log(2)))+1;//log base 2		
		if(pwr%2 == 1){pwr+=1;}
		pxlAraSize = (int) Math.pow(2,pwr);
		dispImgDim = (int)Math.pow(pxlAraSize,.5);
		pxlArray = new int[pxlAraSize];
		maxDispVals=0; minDispVals = 0;
		max2ByteVals=-1110; min2ByteVals = 1110;
		paddedZero = new ArrayList<Integer> ();
		for(int i = origPxlSize; i < pxlAraSize; ++i){
			paddedZero.add(0);
		}
		finishedCalc = false;
		dispImage = pa.createImage(dispImgDim,dispImgDim,PConstants.RGB);
	}

	public void fftTrans(double[] pxls, int w, int h){	
		int scribeIdx = 2;
		//pa.scribeHeader("at start max val : "+Doubles.max(pxls)+"| minVal : "+Doubles.min(pxls),scribeIdx++);
		finishedCalc = false;
		fftRes = fft.transform(pxls, TransformType.FORWARD);
		//dispImage.resize(w, h);
		dispImage.loadPixels();
		int imgAraLen = fftRes.length;
		int halfway = (int)(imgAraLen * .5f)+ (int)(dispImgDim*.5f);
		int tmpIdx = -1, oldIdx, 
				row3Off = (dispImgDim * (dispImgDim/2)),
				colOff2 = (dispImgDim/2),
				newIdx = -1, 
				val;
		for (int row =0; row< dispImgDim; ++row){
			for(int col = 0; col < dispImgDim; ++col){	
				tmpIdx++;
				
				int colPos = (colOff2 + (tmpIdx%dispImgDim)) % dispImgDim;

				oldIdx = tmpIdx;
				newIdx = ((row3Off + (((tmpIdx+colOff2)/dispImgDim)*dispImgDim)) + colPos)%imgAraLen;

				double fftVal =  Math.log(fftRes[oldIdx].abs());
				val = 0xFF & (int)(10*fftVal);
				max2ByteVals = PApplet.max(val, max2ByteVals);
				min2ByteVals  = PApplet.min(val, min2ByteVals);

				//val = 0xFF & (int)(2*mag);
				dispImage.pixels[newIdx] =  0xFF000000 | (0x00FFFFFF & (val<<8 + val<<16 + val));
			}
		}
		
		maxDispVals=PApplet.max(dispImage.pixels); minDispVals = PApplet.min(dispImage.pixels);
		//pa.scribeHeader("max val : "+pa.max(fftResDebug)+"| minVal : "+pa.min(fftResDebug),scribeIdx++);
		pa.scribeHeader("disp pxls max val : "+PApplet.max(dispImage.pixels)+"| minVal : "+PApplet.min(dispImage.pixels),scribeIdx++);
		dispImage.updatePixels();
	    finishedCalc = true;
	}//fftTrans	
	
	public void drawFFT(int xLoc, int yLoc, int width, int height, PImage dispImage){
		pa.scribeHeader("xLoc : "+xLoc+"|YLoc: "+yLoc+" width : "+ (width-xLoc) + " height : " +(height-yLoc), 7);
		pa.pushMatrix();pa.pushStyle();
		pa.beginShape();
		pa.texture(dispImage);
		pa.vertex(xLoc, yLoc, 0, 0);
		pa.vertex(width, yLoc, dispImage.width, 0);
		pa.vertex(width, height, dispImage.width, dispImage.height);
		pa.vertex(xLoc, height, 0, dispImage.height);
		pa.endShape();
			
		pa.popStyle();pa.popMatrix();	
	}
	
	@Override
	public void run() {
		pa.screenShot = pa.get( pa.fftx, pa.ffty,pa.fftw,pa.ffth);
		pa.screenShot.loadPixels();	
		double[] pxls = new double[pxlAraSize];
		//int diffIdx = (fft.pxlAraSize - screenShot.pixels.length)/2;		//center screenshot pxls
//		for(int i =0; i<screenShot.pixels.length; ++i){
//			pxls[i] = (screenShot.pixels[i] >> 16 & 0xFFFF)/65536.0f;
//		}
	
		for(int i =0; i<pa.screenShot.pixels.length; ++i){
			//pxls[diffIdx++] =-1* (screenShot.pixels[i] & 0x00FFFFFF);///256.0f;
			pxls[i] =-1* (pa.screenShot.pixels[i] & 0x00FFFFFF);///256.0f;
		}
		//scribeHeader("max val : "+max(screenShot.pixels)+"| minVal : "+min(screenShot.pixels),6);
		fftTrans(pxls, pa.fftw, pa.ffth);
		// TODO Auto-generated method stub

	}

}
