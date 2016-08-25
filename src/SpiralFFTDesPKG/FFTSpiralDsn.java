package SpiralFFTDesPKG;

import processing.core.*; 

import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.awt.event.KeyEvent;

public class FFTSpiralDsn extends PApplet {
//	
	//**************************** global variables ****************************
	public float sIncr = .01f;	//sIncr is 1/# of frames we wish to display in trajectory display
	public float sW = 1.5f;		//stroke weight
	public Boolean animate=false, 
			drawEdges = false,				//draw edges from each point in spiral
			drawCarrierCurve = false,
					useFFT = true,
			drawPatternCurve = true,
			flipMappedTraj = false,			//flip the direction of the mapped trajectory/spiral on the main spiral (outside or inside)
			endPointsDirty = true,
			shiftKeyPressed = false,
			useDrawnTraj = false,			//use drawn trajectory to determine path 
			useDrawnVels = false;
	
	public ExecutorService th_exec;
	public int numThreadsAvail;	

	//JT drawn trajectory, using drawing speed to infer velocities
	public Boolean drawing = false;		//for drawn trajectory
	public myDrawnObject drawnTraj;						//the 3 objects the user can draw - the rot axis, the left and the right strokes	
	
	public pt[] mainSpiral, spiralEndPts;
	
	public final int fRate = 60;		//frame rate
	
	public PImage screenShot;			//entire screen shot
	public myFFTCalcAndVis fft;			//do real-time fft on image
	
	
	public final int drawnTrajEditWidth = 10;			//width in cntl points of the amount of the drawn trajectory deformed by dragging
	public final float wScale = fRate/5.0f,			//velocity drag scaling	
				trajDragScaleAmt = 100f;					//amt of displacement when dragging drawn trajectory to edit
	public float msClkPtRad = 10;	//radius within which a mouse click will register on a point
	
	public pt spTmplMP;			//midpoint location of spiral template region used to design bounding spiral
	
	public float spTmplOffset, minTmplOff, maxTmplOff;		//offset between end points in template draw region, min and max 
	public int msePickPoint, drawnTrajPickedIdx;	
	public pt spBindA, spBindB,
		mainSpA, mainSpB, mainSpC, mainSpCtr,
		scaleDispPoint;					//primary points for spiral construction
	public float alpha, scaleM,
		minScale = -1.9999f,//.0001f,
		maxScale = 0.9932333f,
		mseSens = 300.0f, 
		templateZoneY;
	public int fftx, ffty,fftw,ffth;		//dims for region to calc fft on
	public float eps = .01f;			//for floating point calcs
	
	//Array that stores all the path Points, once they are scaled
	public pt[] pathBetweenEdges;
	
	
	//**************************** initialization ****************************
	public void settings() {  size(1280, 960, P2D);  smooth(); }				//NEW IN PROCESSING 3.0
	public void setup() {               // executed once at the beginning 		
		//Initialize the singleton
		//instance = this;
		colorMode(RGB, 255);
		drawnTrajPickedIdx = -1;
		msePickPoint = -1;
		spTmplOffset = 400;			//vary this based on scale of drawn arc, to bring the end points together and keep arcs within gray region
		minTmplOff = 10; 		
		maxTmplOff = 400;
		templateZoneY = .7f * this.height;			//region below which the template to modify the spiral binding will be
		spTmplMP = P(.5f * this.width, .85f * this.height);
		mainSpA = P(300,200); 
		mainSpB = P(this.width-300,200); 
		mainSpC = P(this.width * .5f, this.height * .5f); 
		mainSpCtr = spiralCenter(mainSpA, mainSpB,mainSpB, mainSpC);
		scaleM = .01f;
		endPointsDirty = true;

		useFFT = true;
		//point array holding the main spiral's points for display and endpoints for matching pattern - use this so that points are not re-calculated every frame, only when verts A,B,C are moved
		numThreadsAvail = Runtime.getRuntime().availableProcessors();
		System.out.println("# threads : "+ numThreadsAvail);
		th_exec = Executors.newCachedThreadPool();

		mainSpiral = new pt[0];
		spiralEndPts = new pt[0];
		
		updateDrawnSpiral();
		fftx = (int)(max(0,.5f*(int)(this.width - templateZoneY)));
		ffty = 0;
		fftw = 512;//(int)(templateZoneY);//(int)(min(this.width,.5f*(int)(this.width + templateZoneY)));
		ffth = 512;//(int)(templateZoneY);
		screenShot = get( fftx, ffty,fftw,ffth);
		screenShot.loadPixels();
		fft = new myFFTCalcAndVis(this, screenShot);
		frameRate(fRate);             // render 30 frames per second
		
		//animStep = .3f/fRate;			//gives .01 for 30 frames per sec
		pathBetweenEdges = new pt[0];
	}


	//**************************** display current frame ****************************
	public void draw() {      // executed at each frame
		background(white); // clear screen and paints white background		
		
		if(snapPic) beginRecord(PDF,PicturesOutputPath+"/P"+nf(pictureCounter++,3)+".pdf"); 
		//if(animating) {t+=animStep; if(t>=1) {t=1; animating=false;}} 
		drawLegend();
		//array of end points to apply each pattern
		showKeyPt(mainSpA,"A",msClkPtRad);
		showKeyPt(mainSpB,"B",msClkPtRad);
		showKeyPt(mainSpC,"C", msClkPtRad);
		showKeyPt(mainSpCtr,"Ctr",msClkPtRad);
		
		if(endPointsDirty){ buildMainSprial();}
		pushStyle();
		drawMainSpiral();
		popStyle();
//		handle drawn trajectory to give velocity and direction
		pushMatrix();				pushStyle();			
			fill(200,200,200,255);
			rect(0,templateZoneY,this.width,this.height);
			showKeyPt(spBindA,"A1",4);
			showKeyPt(spBindB,"B1",4);		
			if(useDrawnTraj){
				showKeyPt(P(20,templateZoneY+10)," Below is Drawn Pattern",2);
				if(pathBetweenEdges.length > 0) {show(pathBetweenEdges);
				if(drawPatternCurve) {
					pushStyle();drawMainSpiralTraj(((myVariStroke)drawnTraj).getDrawnPtAra(useDrawnVels));popStyle();
				}}
				if(drawnTraj != null){
					((myVariStroke)drawnTraj).drawMe(useDrawnVels);	  			
				}  
			} else {
				showKeyPt(P(20,templateZoneY+10)," Click-Drag in the gray area to modify spiral",2);
				pushStyle();
				drawEmbossedCurve(spBindA, spBindB, scaleM, (flipMappedTraj)? 1: -1);
				popStyle();
				showKeyPt(scaleDispPoint,"Scale is : "+scaleM,2);
				if(drawPatternCurve) {drawMainSpiralSpiral();}
			}
		popStyle();popMatrix();
		if(useFFT){
			if (frameCount % 4 == 1) {
				//th_exec.execute(fft);
				screenShot = get( fftx, ffty,fftw,ffth);
				screenShot.loadPixels();	
				double[] pxls = new double[fft.pxlAraSize];
				//int diffIdx = (fft.pxlAraSize - screenShot.pixels.length)/2;		//center screenshot pxls
//				for(int i =0; i<screenShot.pixels.length; ++i){
//					pxls[i] = (screenShot.pixels[i] >> 16 & 0xFFFF)/65536.0f;
//				}
			
				for(int i =0; i<screenShot.pixels.length; ++i){
					//pxls[diffIdx++] =-1* (screenShot.pixels[i] & 0x00FFFFFF);///256.0f;
					pxls[i] =-1* (screenShot.pixels[i] & 0x00FFFFFF);///256.0f;
				}
				//scribeHeader("max val : "+max(screenShot.pixels)+"| minVal : "+min(screenShot.pixels),6);
				fft.fftTrans(pxls, fftw, ffth);
			}
//			scribeHeader("disp pxls max val : "+fft.maxDispVals+"| minVal : "+fft.minDispVals,4);
//			scribeHeader("calc val pxls max val : "+fft.max2ByteVals+"| minVal : "+fft.min2ByteVals,5);
//			scribeHeader("disp pxls max val : "+PApplet.max(fft.dispImage.pixels)+"| minVal : "+PApplet.min(fft.dispImage.pixels),3);
			fft.drawFFT(0, (int)(this.height*.3f),  (int)(this.width*.3f), (int)(templateZoneY),fft.dispImage);
		}
		if(snapPic) {endRecord(); snapPic=false;} // end saving a .pdf of the screen
		fill(black); displayHeader();
		if(scribeText && !filming) displayFooter(); // shows title, menu, and my face & name 
		if(filming && (animating || change)) saveFrame("FRAMES/F"+nf(frameCounter++,4)+".tif"); // saves a movie frame 
	}  // end of draw()

	public void drawLegend(){
		int idx = 0;		
		fill(useDrawnTraj ? blue : grey); scribeHeader((idx+1)+"- User-Drawn Trajectory - Shift-Click-Drag the Mouse to draw a desired trajectory.  Use drawn velocities for displacement (Press v to toggle) :  " + useDrawnVels,idx+1);idx++; 		  
	}
	
	//draw the curve that begins at one point and ends at the other and consists of the concentric spirals
	//dir is 1 or -1
	public void drawEmbossedCurve (pt begin, pt end, float m, float dir)	{
		pt Mp = P(begin,end); // M is the midpoint of AB - and center of mutual spiral
		drawSpiralOnFlatLine(begin,Mp,m, dir);
		drawSpiralOnFlatLine(end,Mp,m,dir);
	}//drawSpiral
	
	//draw 1 arm of concentric spiral design used to build pattern
	public void drawSpiralOnFlatLine(pt A, pt Mp, float m, float dir){
		float a = PConstants.PI;
		pt B = deriveBFromSpiralCtr(a, m,A,Mp );
		float t = m * d(A,B);
		pt C = MoveByDistanceTowards(B, t, A);
		
		pt curPt, lastPt = A;	
		//1 time calc
		
		float s =spiralScale(A,B,B,C);
		pt G = spiralCenter( dir*a, s, A, B); 
		float maxIter = PApplet.abs(7/(1-s));
		maxIter = (maxIter > 400 ? 400 : maxIter);
		//System.out.println("scale in drawspiral "+ s + "  |  " +  maxIter);
		for (float i=-1; i<maxIter; i+=sIncr) {			
			curPt = L(G,R(B,dir*i*a,G),PApplet.pow(s,i));
			edge(curPt,lastPt);
			lastPt = curPt;
		}
	}//drawSpiral
		
	public void showKeyPt(pt a, String s, float rad){
		a.show(rad);
		a.label(s);
	}
	//build main spiral every time cntl point positions are changed
	public void buildMainSprial(){
		//mainSpA, mainSpB, mainSpC, mainSpCtr
		float lclalpha =spiralAngle(mainSpA,mainSpB,mainSpB,mainSpC); 
		float lclScaleM =spiralScale(mainSpA,mainSpB,mainSpB,mainSpC);	
		float maxIter = PApplet.abs(7/(1-lclScaleM));
		maxIter = (maxIter > 200 ? 200 : maxIter);		//20k pts - need to make this smaller, most of these points are wasted
		ArrayList<pt> msp = new ArrayList<pt>();
		ArrayList<pt> sep = new ArrayList<pt>();
		mainSpiral = new pt[round(((maxIter+1)/sIncr) +1)];
		spiralEndPts = new pt[round(mainSpiral.length*sIncr)+1];
		pt c;
		for(float s=-1; s<=maxIter; s+=sIncr) {				
			c = L(mainSpCtr,R(mainSpB,s*lclalpha,mainSpCtr),pow(lclScaleM,s));
			msp.add(c);
			if(abs(s-((int)s)) < eps){	sep.add(c);	}						
		}
		mainSpiral = msp.toArray(new pt[0]);
		spiralEndPts = sep.toArray(new pt[0]);
		endPointsDirty = false;
	}
	//draw the precomputed points in the array
	public void drawMainSpiral(){
		if(drawCarrierCurve){
			stroke(red);strokeWeight(sW);
			curve(mainSpiral);
		}
		if(drawEdges){
			stroke(green);strokeWeight(sW);
			showNoClose(spiralEndPts);
		}
	}
	
	//draw spiral embossed on main spiral
	public void drawMainSpiralSpiral(){		
		int numCurves = spiralEndPts.length;
		stroke(purple);strokeWeight(sW);
		float dir;
		if(flipMappedTraj){		dir= 1; } else {dir=-1;}
		for(int i =0; i< numCurves-1; ++i){			
			drawEmbossedCurve(spiralEndPts[i], spiralEndPts[i+1], scaleM, dir);				
		}
	}//drawMainSpiralSpiral
	//draw drawn trajectory embossed on spiral
	public void drawMainSpiralTraj(pt[] src){	
		pt[] res = new pt[src.length];
		int numCurves = spiralEndPts.length;
		stroke(blue);strokeWeight(sW);
		int a,b;
		if(flipMappedTraj){		a = 0; b= 1; } else {a = 1; b= 0;}
		for(int i =0; i< numCurves-1; ++i){			
			res =buildDupeCurveBetween(src,spiralEndPts[i+a], spiralEndPts[i+b]);  
			curve(res);
		}	
		//curve(res);
	}//drawMainSpiralTraj
	//this function derives the b point of a spiral tri-point system, given the angle, scale, A point and ctr
	public pt deriveBFromSpiralCtr(float angle, float scale, pt A, pt Ctr){
		  float c=PApplet.cos(angle), s=PApplet.sin(angle);
		  float D = PApplet.sq((c*scale)-1)+PApplet.sq(s*scale);
		  //precalc constants for B derivation
		  float T1 = c*scale*A.x - s*scale*A.y,
				T2 = c*scale*A.y + s*scale*A.x;
		  float T3 = (c*scale-1), T4 = s*scale;
		  float dFx = Ctr.x * D, dFy = Ctr.y * D;
		  float T3sqMnT4sq = (T3*T3)-(T4*T4),			//precalc
	      By = ((T3*T3*T2) - (T2*T4*T4) + (T4 * dFx) - (T3 * dFy))/T3sqMnT4sq,
		  Bx =  T1 + (T4*T2/T3) - (By*T4/T3) - (dFx/T3);
		  
		  return P(Bx,By);
	}

	//**************************** user actions ****************************
	public void keyPressed() { // executed each time a key is pressed: sets the "keyPressed" and "key" state variables, 
	                    // till it is released or another key is pressed or released
		if(key=='?') scribeText=!scribeText; // toggle display of help text and authors picture
		if(key=='!') snapPicture(); // make a picture of the canvas and saves as .jpg image
		if(key=='`') snapPic=true; // to snap an image of the canvas and save as zoomable a PDF
		if(key=='~') { filming=!filming; } // filming on/off capture frames into folder FRAMES 
		if(key=='1') {useDrawnTraj = !useDrawnTraj; if(useDrawnTraj){useDrawnVels=false; }	}			//use trajectory drawn by user to determine path and velocity
		if(key=='2') {drawEdges = !drawEdges;}								//draw edges from each point in spiral
		if(key=='3') {drawCarrierCurve = !drawCarrierCurve;}				//draw edges from each point in spiral
		if(key=='4') {drawPatternCurve = !drawPatternCurve;}				//draw edges from each point in spiral
		if(key=='5') {flipMappedTraj = !flipMappedTraj;}					//flip whether the drawn trajectory is on the inside or the outside of the base spiral
		if(key=='v') {useDrawnVels = !useDrawnVels; rebuildDrawnTraj();}			//use velocities of drawn trajectory for movement along with angle change
		if(key=='f') {useFFT = !useFFT;}			//display fft info
		
		if(key=='Q') exit();  // quit application
		change=true; // to make sure that we save a movie frame each time something changes
		if((!shiftKeyPressed)&&(key==CODED)){shiftKeyPressed = (keyCode  == KeyEvent.VK_SHIFT);}
	}
	public void keyReleased(){if((shiftKeyPressed)&&(key==CODED)){ if(keyCode == KeyEvent.VK_SHIFT){shiftKeyPressed = false;if (drawing) {endDrawObj();}}}}	
	
	public void mousePressed() {  // executed when the mouse is pressed
		if((useDrawnTraj) && (shiftKeyPressed)){
			//begin drawing object
			startDrawObj(); 
		}
		else{
			//Pick point to drag
			pt mse = Mouse();
			if((useDrawnTraj) && (mse.y > templateZoneY)){	//drawing in template zone - editing curve
				System.out.println("click in click zone");
				float[] distToPt = new float[1];			//using array as pointer, passing by reference
				pt[] pts = ((myVariStroke)drawnTraj).getDrawnPtAra(useDrawnVels);
				int pIdx = drawnTraj.findClosestPt(mse, distToPt,pts);
				if(distToPt[0] < msClkPtRad){//close enough to mod
					System.out.println("dist : "+distToPt[0] + " idx : "+ pIdx);
					drawnTrajPickedIdx = pIdx;					
				}
			} else{	
				//ugh, should put these points in an array
				if(d(mse,mainSpA) < msClkPtRad){msePickPoint = 0;}
				else if(d(mse,mainSpB) < msClkPtRad){msePickPoint = 1;}
				else if(d(mse,mainSpC) < msClkPtRad){msePickPoint = 2;}
				else if(d(mse,mainSpCtr) < msClkPtRad){msePickPoint = 3;}
			}
		}
		change=true;
	 }

	public void mouseDragged() {
		pt mse = Mouse();
		if ((drawing) && (shiftKeyPressed)){
			drawnTraj.addPt(Mouse());
		} else if((useDrawnTraj) && (drawnTrajPickedIdx != -1)){	
			//needs to be before templateZoneY check
			((myVariStroke)drawnTraj).handleMouseDrag( mse,V((mouseX-pmouseX)/mseSens,(mouseY-pmouseY)/mseSens),drawnTrajPickedIdx);
		} else if (mse.y > templateZoneY){//modify scale of ornament here, or modify drawn trajectory	
			msePickPoint = -1;			
			scaleM -=(mouseY-pmouseY)/mseSens; 
			scaleM +=(mouseX-pmouseX)/mseSens;
			scaleM = ((scaleM < maxScale) && (scaleM > minScale) ? scaleM  : (scaleM < minScale ?  minScale : maxScale));
			updateDrawnSpiral();
		} else {//modify main points
			//ugh, should put these points in an array
			if(msePickPoint != -1){
				pt del = P((mouseX-pmouseX),(mouseY-pmouseY));
				if(msePickPoint == 3){//center, move all points by delta
					mainSpA.add(del);
					mainSpB.add(del);
					mainSpC.add(del);	
					mainSpCtr.add(del);
					endPointsDirty = true;
				}
				else if(msePickPoint == 0){mainSpA.add(del);updateMainSpiral();}
				else if(msePickPoint == 1){mainSpB.add(del);updateMainSpiral();}
				else if(msePickPoint == 2){
					mainSpC.add(del);
					float distAB = d(mainSpA,mainSpB);
					if(d(mainSpC,mainSpB) > distAB){mainSpC = P(mainSpB,distAB,U(mainSpB,mainSpC));}//cap c so it can never be further from b than a
					updateMainSpiral();
				}
			}
		}
		change=true;
	}//  mouseDragged

	public void mouseReleased(){
		if((msePickPoint == 0) || (msePickPoint == 1) || (msePickPoint == 2)){	updateMainSpiral();}
		msePickPoint = -1;
		if(drawnTrajPickedIdx != -1){
			drawnTrajPickedIdx = -1;
			drawnTraj.remakeDrawnTraj(useDrawnVels);
			rebuildDrawnTraj();			
		}
		if ((drawing) && (shiftKeyPressed)){endDrawObj(); return;}	
	}	
	
	public void startDrawObj(){
		drawing = true;
		drawnTraj = new myVariStroke(this);
		drawnTraj.startDrawing();
	}
	
	public void rebuildDrawnTraj(){
		//Once edge is drawn
		if(drawnTraj != null){
			if(drawnTraj.flags[drawnTraj.isMade]){				  
				//Initialize the array that stores the path
				pt a,b;
				if(flipMappedTraj){		a = spBindA; b= spBindB;} else {a = spBindB; b= spBindA;}
				if(useDrawnVels){pathBetweenEdges =((myVariStroke)drawnTraj).moveVelCurveToEndPoints(a, b); }
				else{pathBetweenEdges =((myVariStroke)drawnTraj).moveCntlCurveToEndPoints(a, b);  	}
			}
		}		
	}
	
	public void endDrawObj(){
		drawnTraj.addPt(Mouse());
		drawnTraj.finalizeDrawing(true);	
		
		rebuildDrawnTraj();
		drawing = false;
		//applyTrajTo
	}
	
	public void updateDrawnSpiral(){
		if(scaleM < 0){spTmplOffset = minTmplOff + .5f*(2+scaleM) * (maxTmplOff - minTmplOff);} else {spTmplOffset = minTmplOff + (1-scaleM) * (maxTmplOff - minTmplOff);}
		spBindA = P(spTmplMP.x - spTmplOffset,spTmplMP.y);
		spBindB = P(spTmplMP.x + spTmplOffset,spTmplMP.y);
		scaleDispPoint = P(spTmplMP.x + 200,spTmplMP.y * .95f);	
		if(useDrawnTraj){
			rebuildDrawnTraj();
		}
	}
	public void updateMainSpiral(){
		//if any point has been moved, modify spiral and center accordingly
		  float a =spiralAngle(mainSpA,mainSpB,mainSpB,mainSpC); 
		  float s =spiralScale(mainSpA,mainSpB,mainSpB,mainSpC);
		  mainSpCtr = spiralCenter(a, s, mainSpA, mainSpB); 
		  endPointsDirty = true;
	}
	
	
	public pt[] buildDupeCurveBetween(pt[] srcCurve, pt startPt, pt endPt){
		pt[] destCurve = new pt[srcCurve.length];
		//drawn curve params
		if(srcCurve.length == 0){return destCurve;}
		pt origin = srcCurve[0];
		pt end = srcCurve[srcCurve.length - 1];
		
		//edge params		
		vec drawnAxis = V(origin, end);
		vec edgeAxis =  V(startPt, endPt);		//angle between these two is the angle to rotate everyone
		
		//transformation params
		vec dispToStart = V(origin, startPt);			//displacement vector between start of drawn curve and edge 1.
		float alpha = angle(drawnAxis,edgeAxis);			//angle to rotate everyone
		float scaleRatio = n(edgeAxis)/n(drawnAxis);	//ratio of distance from start to finish of drawn traj to distance between edges - multiply all elements in drawn traj by this
		
		//displace to align with start
		for(int ptItr = 0; ptItr < srcCurve.length ; ++ptItr){
			destCurve[ptItr] = P(srcCurve[ptItr], dispToStart);
		}
		//displace every point to be scaled distance from start of curve equivalent to scale of edge distances to drawn curve
		for(int ptItr = 1; ptItr < srcCurve.length ; ++ptItr){
			destCurve[ptItr] = P(destCurve[0],scaleRatio, V(destCurve[0],destCurve[ptItr]));//start point displaced by scaleRatio * vector from start to original location of pt
		}
		//rotate every point around destCurve[0] by alpha
		for(int ptItr = 1; ptItr < srcCurve.length ; ++ptItr){
			destCurve[ptItr] = R(destCurve[ptItr],alpha,destCurve[0]);
		}
		return destCurve;
		
	}//fitCurveBetweenAlt
	
	//**************************** text for name, title and help  ****************************
	String title ="CA 2015 P2: Pattern Design", 
	       name ="John Turner",
	       menu="?:(show/hide) help, f: fft freq window, Q:quit",
	       guide="1: drawn/spiral embossing, 2: draw edges, 3: draw carrier curve,  4 : draw embossed pattern, 5 : flip direction of drawn curve, v : drawn pos/drawn vels used for embossed curve"; // help info
	
	public void drawObject(pt P, vec V) {
	  beginShape(); 
	    v(P(P(P,1,V),0.25f,R(V)));
	    v(P(P(P,1,V),-0.25f,R(V)));
	    v(P(P(P,-1,V),-0.25f,R(V)));
	    v(P(P(P,-1,V),0.25f,R(V))); 
	  endShape(CLOSE);
	  }
	  
	public float timeWarp(float f) {return sq(sin(f*PI/2));}
	//************************************************************************
	//**** CIRCLES
	//************************************************************************
	// create 
	public float circumRadius (pt A, pt B, pt C) {float a=d(B,C), b=d(C,A), c=d(A,B), s=(a+b+c)/2, d=sqrt(s*(s-a)*(s-b)*(s-c)); return a*b*c/4/d;} // radiusCircum(A,B,C): radius of circumcenter
	public pt CircumCenter (pt A, pt B, pt C) {vec AB = V(A,B); vec AC = R(V(A,C)); 
	   return P(A,1.f/2/dot(AB,AC),W(-n2(AC),R(AB),n2(AB),AC)); }; // CircumCenter(A,B,C): center of circumscribing circle, where medians meet)
	
	// display 
	public void drawCircle(int n) {  
	  float x=1, y=0; float a=TWO_PI/n, t=tan(a/2), s=sin(a); 
	  beginShape(); for (int i=0; i<n; i++) {x-=y*t; y+=x*s; x-=y*t; vertex(x,y);} endShape(CLOSE);}
	
	
	public void showArcThrough (pt A, pt B, pt C) {
	  if (abs(dot(V(A,B),R(V(A,C))))<0.01f*d2(A,C)) {edge(A,C); return;}
	   pt O = CircumCenter ( A,  B,  C); 
	   float r=d(O,A);
	   vec OA=V(O,A), OB=V(O,B), OC=V(O,C);
	   float b = angle(OA,OB), c = angle(OA,OC); 
	   if(0<c && c<b || b<0 && 0<c)  c-=TWO_PI; 
	   else if(b<c && c<0 || c<0 && 0<b)  c+=TWO_PI; 
	   beginShape(); v(A); for (float t=0; t<1; t+=0.01f) v(R(A,t*c,O)); v(C); endShape();
	   }
	
	public pt pointOnArcThrough (pt A, pt B, pt C, float t) { // July 2011
	  if (abs(dot(V(A,B),R(V(A,C))))<0.001f*d2(A,C)) {edge(A,C); return L(A,C,t);}
	   pt O = CircumCenter ( A,  B,  C); 
	   float r=(d(O,A) + d(O,B)+ d(O,C))/3;
	   vec OA=V(O,A), OB=V(O,B), OC=V(O,C);
	   float b = angle(OA,OB), c = angle(OA,OC); 
	   if(0<b && b<c) {}
	   else if(0<c && c<b) {b=b-TWO_PI; c=c-TWO_PI;}
	   else if(b<0 && 0<c) {c=c-TWO_PI;}
	   else if(b<c && c<0) {b=TWO_PI+b; c=TWO_PI+c;}
	   else if(c<0 && 0<b) {c=TWO_PI+c;}
	   else if(c<b && b<0) {}
	   return R(A,t*c,O);
	   }

	//*****************************************************************************
	// TITLE:         GEOMETRY UTILITIES IN 2D  
	// DESCRIPTION:   Classes and functions for manipulating points, vectors, edges, triangles, quads, frames, and circular arcs  
	// AUTHOR:        Prof Jarek Rossignac
	// DATE CREATED:  September 2009
	// EDITS:         Revised July 2011
	//*****************************************************************************
	//************************************************************************
	//**** POINT CLASS
	//************************************************************************
	public class pt { float x=0,y=0; 
	  // CREATE
	  pt () {}
	  pt (float px, float py) {x = px; y = py;};
	
	  // MODIFY
	  public pt setTo(float px, float py) {x = px; y = py; return this;};  
	  public pt set(pt P) {x = P.x; y = P.y; return this;}; 
	  public pt setTo(pt P) {x = P.x; y = P.y; return this;}; 
	  public pt setToMouse() { x = mouseX; y = mouseY;  return this;}; 
	  public pt add(float u, float v) {x += u; y += v; return this;}                       // P.add(u,v): P+=<u,v>
	  public pt add(pt P) {x += P.x; y += P.y; return this;};                              // incorrect notation, but useful for computing weighted averages
	  public pt add(float s, pt P)   {x += s*P.x; y += s*P.y; return this;};               // adds s*P
	  public pt add(vec V) {x += V.x; y += V.y; return this;}                              // P.add(V): P+=V
	  public pt add(float s, vec V) {x += s*V.x; y += s*V.y; return this;}                 // P.add(s,V): P+=sV
	  public pt translateTowards(float s, pt P) {x+=s*(P.x-x);  y+=s*(P.y-y);  return this;};  // transalte by ratio s towards P
	  public pt scale(float u, float v) {x*=u; y*=v; return this;};
	  public pt scale(float s) {x*=s; y*=s; return this;}                                  // P.scale(s): P*=s
	  public pt scale(float s, pt C) {x*=C.x+s*(x-C.x); y*=C.y+s*(y-C.y); return this;}    // P.scale(s,C): scales wrt C: P=L(C,P,s);
	  public pt rotate(float a) {float dx=x, dy=y, c=cos(a), s=sin(a); x=c*dx+s*dy; y=-s*dx+c*dy; return this;};     // P.rotate(a): rotate P around origin by angle a in radians
	  public pt rotate(float a, pt G) {float dx=x-G.x, dy=y-G.y, c=cos(a), s=sin(a); x=G.x+c*dx+s*dy; y=G.y-s*dx+c*dy; return this;};   // P.rotate(a,G): rotate P around G by angle a in radians
	  public pt rotate(float s, float t, pt G) {float dx=x-G.x, dy=y-G.y; dx-=dy*t; dy+=dx*s; dx-=dy*t; x=G.x+dx; y=G.y+dy;  return this;};   // fast rotate s=sin(a); t=tan(a/2); 
	  public pt moveWithMouse() { x += mouseX-pmouseX; y += mouseY-pmouseY;  return this;}; 
	     
	  // DRAW , WRITE
	  public pt write() {print("("+x+","+y+")"); return this;};  // writes point coordinates in text window
	  public pt v() {vertex(x,y); return this;};  // used for drawing polygons between beginShape(); and endShape();
	  public pt show(float r) {ellipse(x, y, 2*r, 2*r); return this;}; // shows point as disk of radius r
	  public void v(pt P) {vertex(P.x,P.y);};                                           // vertex for shading or drawing
	  public pt show() {show(3); return this;}; // shows point as small dot
	  public pt label(String s, float u, float v) {fill(black); text(s, x+u, y+v); noFill(); return this; };
	  public pt label(String s, vec V) {fill(black); text(s, x+V.x, y+V.y); noFill(); return this; };
	  public pt label(String s) {label(s,5,4); return this; };
	  public String toString(){return "{x: "+x+", y: "+y+"}";};
	  } // end of pt class
	
	//************************************************************************
	//**** VECTOR CLASS
	//************************************************************************
	public class vec { float x=0,y=0; 
	 // CREATE
	  vec () {};
	  vec (float px, float py) {x = px; y = py;};
	 
	 // MODIFY
	  public vec setTo(float px, float py) {x = px; y = py; return this;}; 
	  public vec setTo(vec V) {x = V.x; y = V.y; return this;}; 
	  public vec zero() {x=0; y=0; return this;}
	  public vec scaleBy(float u, float v) {x*=u; y*=v; return this;};
	  public vec scaleBy(float f) {x*=f; y*=f; return this;};
	  public vec mul(float f) {x*=f; y*=f; return this;}	  
	  public vec reverse() {x=-x; y=-y; return this;};
	  public vec divideBy(float f) {x/=f; y/=f; return this;};
	  public vec normalize() {float n=sqrt(sq(x)+sq(y)); if (n>0.000001f) {x/=n; y/=n;}; return this;};
	  public vec add(float u, float v) {x += u; y += v; return this;};
	  public vec add(vec V) {x += V.x; y += V.y; return this;};   
	  public vec add(float s, vec V) {x += s*V.x; y += s*V.y; return this;};   
	  public vec rotateBy(float a) {float xx=x, yy=y; x=xx*cos(a)-yy*sin(a); y=xx*sin(a)+yy*cos(a); return this;};
	  public vec left() {float m=x; x=-y; y=m; return this;};
	 
	  // OUTPUT VEC
	  public vec clone() {return(new vec(x,y));}; 
	  public String toString(){return "{x: "+x+", y: "+y+"}";};
	
	  // OUTPUT TEST MEASURE
	  public float norm() {return(sqrt(sq(x)+sq(y)));}
	  public boolean isNull() {return((abs(x)+abs(y)<0.000001f));}
	  public float angle() {return(atan2(y,x)); }
	
	  // DRAW, PRINT
	  public void write() {println("<"+x+","+y+">");};
	  public void showAt (pt P) {line(P.x,P.y,P.x+x,P.y+y); }; 
	  public void showArrowAt (pt P) {line(P.x,P.y,P.x+x,P.y+y); 
	      float n=min(this.norm()/10.f,height/50.f); 
	      pt Q=P(P,this); 
	      vec U = S(-n,U(this));
	      vec W = S(.3f,R(U)); 
	      beginShape(); Q.add(U).add(W).v(); Q.v(); Q.add(U).add(M(W)).v(); endShape(CLOSE); }; 
	  public void label(String s, pt P) {P(P).add(0.5f,this).add(3,R(U(this))).label(s); };
	  } // end vec class
	
	//************************************************************************
	//**** POINTS FUNCTIONS
	//************************************************************************
	// create 
	public pt P() {return P(0,0); };                                                                            // make point (0,0)
	public pt P(float x, float y) {return new pt(x,y); };                                                       // make point (x,y)
	public pt P(pt P) {return P(P.x,P.y); };                                                                    // make copy of point A
	public pt P(pt O, float x, vec I, float y, vec J) {return P(O.x+x*I.x+y*J.x,O.y+x*I.y+y*J.y);}  // O+xI+yJ
	public pt Mouse() {return P(mouseX,mouseY);};                                                                 // returns point at current mouse location
	public pt Pmouse() {return P(pmouseX,pmouseY);};                                                              // returns point at previous mouse location
	public pt ScreenCenter() {return P(width/2,height/2);}                                                        //  point in center of  canvas
	public pt P(vec V) { return P(V.x,V.y);}
	// transform 
	public pt R(pt Q, float a) {float dx=Q.x, dy=Q.y, c=cos(a), s=sin(a); return new pt(c*dx+s*dy,-s*dx+c*dy); };  // Q rotated by angle a around the origin
	public pt R(pt Q, float a, pt C) {float dx=Q.x-C.x, dy=Q.y-C.y, c=cos(a), s=sin(a); return P(C.x+c*dx-s*dy, C.y+s*dx+c*dy); };  // Q rotated by angle a around point C
	public pt P(pt P, vec V) {return P(P.x + V.x, P.y + V.y); }                                                 //  P+V (P transalted by vector V)
	public pt P(pt P, float s, vec V) {return P(P,W(s,V)); }                                                    //  P+sV (P transalted by sV)
	public pt MoveByDistanceTowards(pt P, float d, pt Q) { return P(P,d,U(V(P,Q))); };                          //  P+dU(PQ) (transLAted P by *distance* s towards Q)!!!
	
	// average 
	public pt P(pt A, pt B) {return P((A.x+B.x)/2.0f,(A.y+B.y)/2.0f); };                                          // (A+B)/2 (average)
	public pt P(pt A, pt B, pt C) {return P((A.x+B.x+C.x)/3.0f,(A.y+B.y+C.y)/3.0f); };                            // (A+B+C)/3 (average)
	public pt P(pt A, pt B, pt C, pt D) {return P(P(A,B),P(C,D)); };                                            // (A+B+C+D)/4 (average)
	
	// weighted average 
	public pt P(float a, pt A) {return P(a*A.x,a*A.y);}                                                      // aA  
	public pt P(float a, pt A, float b, pt B) {return P(a*A.x+b*B.x,a*A.y+b*B.y);}                              // aA+bB, (a+b=1) 
	public pt P(float a, pt A, float b, pt B, float c, pt C) {return P(a*A.x+b*B.x+c*C.x,a*A.y+b*B.y+c*C.y);}   // aA+bB+cC 
	public pt P(float a, pt A, float b, pt B, float c, pt C, float d, pt D){return P(a*A.x+b*B.x+c*C.x+d*D.x,a*A.y+b*B.y+c*C.y+d*D.y);} // aA+bB+cC+dD 
	
	// LERP
	public pt L(pt A, pt B, float t) {return P(A.x+t*(B.x-A.x),A.y+t*(B.y-A.y));}
	public pt L(pt A, float t, pt B) {return L(A,B,t);}								//just a different config for imported classes
			    
	// measure 
	public boolean isSame(pt A, pt B) {return (A.x==B.x)&&(A.y==B.y) ;}                                         // A==B
	public boolean isSame(pt A, pt B, float e) {return ((abs(A.x-B.x)<e)&&(abs(A.y-B.y)<e));}                   // ||A-B||<e
	public float d(pt P, pt Q) {return sqrt(d2(P,Q));  };                                                       // ||AB|| (Distance)
	public float d2(pt P, pt Q) {return sq(Q.x-P.x)+sq(Q.y-P.y); };                                             // AB*AB (Distance squared)
	
	//************************************************************************
	//**** VECTOR FUNCTIONS
	//************************************************************************
	// create 
	public vec V(vec V) {return new vec(V.x,V.y); };                                                             // make copy of vector V
	public vec V(pt P) {return new vec(P.x,P.y); };                                                              // make vector from origin to P
	public vec V(float x, float y) {return new vec(x,y); };                                                      // make vector (x,y)
	public vec V(pt P, pt Q) {return new vec(Q.x-P.x,Q.y-P.y);};                                                 // PQ (make vector Q-P from P to Q
	public vec U(vec V) {float n = n(V); if (n==0) return new vec(0,0); else return new vec(V.x/n,V.y/n);};      // V/||V|| (Unit vector : normalized version of V)
	public vec U(pt P, pt Q) {return U(V(P,Q));};                                                                // PQ/||PQ| (Unit vector : from P towards Q)
	public vec MouseDrag() {return new vec(mouseX-pmouseX,mouseY-pmouseY);};                                      // vector representing recent mouse displacement
	
	// weighted sum 
	public vec W(float s,vec V) {return V(s*V.x,s*V.y);}                                                      // sV
	public vec W(vec U, vec V) {return V(U.x+V.x,U.y+V.y);}                                                   // U+V 
	public vec W(vec U,float s,vec V) {return W(U,S(s,V));}                                                   // U+sV
	public vec W(float u, vec U, float v, vec V) {return W(S(u,U),S(v,V));}                                   // uU+vV ( Linear combination)
	
	// transformed 
	public vec R(vec V) {return new vec(-V.y,V.x);};                                                             // V turned right 90 degrees (as seen on screen)
	public vec R(vec V, float a) {float c=cos(a), s=sin(a); return(new vec(V.x*c-V.y*s,V.x*s+V.y*c)); };                                     // V rotated by a radians
	public pt R(pt P, float a, vec I, vec J, pt G) {float x=dot(V(G,P),I), y=dot(V(G,P),J); float c=cos(a), s=sin(a); return P(P,x*c-x-y*s,I,x*s+y*c-y,J); }; 
	public cntlPt R(cntlPt P, float a, vec I, vec J, pt G) {float x=dot(V(G,P),I), y=dot(V(G,P),J); float c=cos(a), s=sin(a); return new cntlPt(this, P(P,x*c-x-y*s,I,x*s+y*c-y,J), P.r, P.w); }; 
	
	public vec S(float s,vec V) {return new vec(s*V.x,s*V.y);};                                                  // sV
	public vec Reflection(vec V, vec N) { return W(V,-2.f*dot(V,N),N);};                                          // reflection
	public vec M(vec V) { return V(-V.x,-V.y); }                                                                  // -V
	
	// Interpolation 
	public vec L(vec U, vec V, float s) {return new vec(U.x+s*(V.x-U.x),U.y+s*(V.y-U.y));};                      // (1-s)U+sV (Linear interpolation between vectors)
	public vec L(vec U, float s, vec V) {return new vec(U.x+s*(V.x-U.x),U.y+s*(V.y-U.y));};                      // alt call for imported functions

	public vec S(vec U, vec V, float s) {float a = angle(U,V); vec W = R(U,s*a); float u = n(U), v=n(V); return W(pow(v/u,s),W); } // steady interpolation from U to V
	public vec S(vec U, float s, vec V) {float a = angle(U,V); vec W = R(U,s*a); float u = n(U), v=n(V); return W(pow(v/u,s),W); } // alt call for imported functions
	
	// measure 
	public float dot(vec U, vec V) {return U.x*V.x+U.y*V.y; }                                                     // dot(U,V): U*V (dot product U*V)
	public float det(vec U, vec V) {return dot(R(U),V); }                                                         // det | U V | = scalar cross UxV 
	public float n(vec V) {return sqrt(dot(V,V));};                                                               // n(V): ||V|| (norm: length of V)
	public float n2(vec V) {return sq(V.x)+sq(V.y);};                                                             // n2(V): V*V (norm squared)
	public boolean parallel (vec U, vec V) {return dot(U,R(V))==0; }; 
	
	public float angle (vec U, vec V) {return atan2(det(U,V),dot(U,V)); };                                   // angle <U,V> (between -PI and PI)
	public float angle(vec V) {return(atan2(V.y,V.x)); };                                                       // angle between <1,0> and V (between -PI and PI)
	public float angle(pt A, pt B, pt C) {return  angle(V(B,A),V(B,C)); }                                       // angle <BA,BC>
	public float turnAngle(pt A, pt B, pt C) {return  angle(V(A,B),V(B,C)); }                                   // angle <AB,BC> (positive when right turn as seen on screen)
	public int toDeg(float a) {return PApplet.parseInt(a*180/PI);}                                                           // convert radians to degrees
	public float toRad(float a) {return(a*PI/180);}                                                             // convert degrees to radians 
	public float positive(float a) { if(a<0) return a+TWO_PI; else return a;}                                   // adds 2PI to make angle positive
	
	// SLERP
	public vec slerp(vec U, float t, vec V) {float a = angle(U,V); float b=sin((1.f-t)*a),c=sin(t*a),d=sin(a); return W(b/d,U,c/d,V); } // UNIT vectors ONLY!
	
	//************************************************************************
	//**** DISPLAY
	//************************************************************************
	// point / polygon
	public void show(pt P, float r) {ellipse(P.x, P.y, 2*r, 2*r);};                                             // draws circle of center r around P
	public void show(pt P) {ellipse(P.x, P.y, 6,6);};                                                           // draws small circle around point
	public void show(pt[] ara) {beginShape(); for(int i=0;i<ara.length;++i){v(ara[i]);} endShape(CLOSE);};                    
	public void showNoClose(pt[] ara) {beginShape(); for(int i=0;i<ara.length;++i){v(ara[i]);} endShape();};                    
	public void curveVertex(pt P) {curveVertex((float)P.x,(float)P.y);};                                           // curveVertex for shading or drawing
	public void curve(pt[] ara) {if(ara.length == 0){return;}beginShape(); curveVertex(ara[0]);for(int i=0;i<ara.length;++i){curveVertex(ara[i]);} curveVertex(ara[ara.length-1]);endShape();};                      // volume of tet 
	
	//draw a circle - JT
	void circle(pt P, float r, vec I, vec J, int n) {pt[] pts = new pt[n];pts[0] = P(P,r,U(I));float a = (2*PI)/(1.0f*n);for(int i=1;i<n;++i){pts[i] = R(pts[i-1],a,J,I,P);}pushMatrix(); pushStyle();noFill(); show(pts);popStyle();popMatrix();}; // render sphere of radius r and center P
	
	// edge / arrow
	public void edge(pt P, pt Q) {line(P.x,P.y,Q.x,Q.y); };                                                      // draws edge (P,Q)
	public void arrow(pt P, pt Q) {arrow(P,V(P,Q)); }                                                            // draws arrow from P to Q
	public void show(pt P, vec V) {line(P.x,P.y,P.x+V.x,P.y+V.y); }                                              // show V as line-segment from P 
	public void show(pt P, float s, vec V) {show(P,S(s,V));}                                                     // show sV as line-segment from P 
	public void arrow(pt P, float s, vec V) {arrow(P,S(s,V));}                                                   // show sV as arrow from P 
	public void arrow(pt P, vec V, String S) {arrow(P,V); P(P(P,0.70f,V),15,R(U(V))).label(S,V(-5,4));}       // show V as arrow from P and print string S on its side
	public void arrow(pt P, vec V) {show(P,V);  float n=n(V); if(n<0.01f) return; float s=max(min(0.2f,20.f/n),6.f/n);       // show V as arrow from P 
	     pt Q=P(P,V); vec U = S(-s,V); vec W = R(S(.3f,U)); beginShape(); v(P(P(Q,U),W)); v(Q); v(P(P(Q,U),-1,W)); endShape(CLOSE);}; 
	
	// triangle, polygon
	public void v(pt P) {vertex(P.x,P.y);};                                                                     // vertex for drawing polygons between beginShape() and endShape()
	public void v(pt P, float u, float v) { vertex(P.x, P.y, u, v); }
	public void show(pt A, pt B, pt C)  {beginShape();  A.v(); B.v(); C.v(); endShape(CLOSE);}                   // render triangle A, B, C
	public void show(pt A, pt B, pt C, pt D)  {beginShape();  A.v(); B.v(); C.v(); D.v(); endShape(CLOSE);}      // render quad A, B, C, D
	
	// text
	public void label(pt P, String S) {text(S, P.x-4,P.y+6.5f); }                                                 // writes string S next to P on the screen ( for example label(P[i],str(i));)
	public void label(pt P, vec V, String S) {text(S, P.x-3.5f+V.x,P.y+7+V.y); }                                  // writes string S at P+V
	//************************************************************************
	//**** SPIRAL
	//************************************************************************
	public pt PtOnSpiral(pt A, pt B, pt C, float t) {
	  float a =spiralAngle(A,B,B,C); 
	  float s =spiralScale(A,B,B,C);
	  pt G = spiralCenter(a, s, A, B); 
	  return L(G,R(B,t*a,G),pow(s,t));
	  }
	
	public pt spiralPt(pt A, pt G, float s, float a) {return L(G,R(A,a,G),s);}  
	public pt spiralPt(pt A, pt G, float s, float a, float t) {return L(G,R(A,t*a,G),pow(s,t));} 
	public pt spiralCenter(pt A, pt B, pt C, pt D) { // computes center of spiral that takes A to C and B to D
	  float a = spiralAngle(A,B,C,D); 
	  float z = spiralScale(A,B,C,D);
	  return spiralCenter(a,z,A,C);
	  }
	public float spiralAngle(pt A, pt B, pt C, pt D) {return angle(V(A,B),V(C,D));}
	public float spiralScale(pt A, pt B, pt C, pt D) {return d(C,D)/d(A,B);}
	public pt spiralCenter(float a, float z, pt A, pt C) {
	  float c=cos(a), s=sin(a);
	  float D = sq(c*z-1)+sq(s*z);
	  float ex = c*z*A.x - C.x - s*z*A.y;
	  float ey = c*z*A.y - C.y + s*z*A.x;
	  float x=(ex*(c*z-1) + ey*s*z) / D;
	  float y=(ey*(c*z-1) - ex*s*z) / D;
	  return P(x,y);
	  }	  
	  
	// ************************************************************************ COLORS 
	int black=0xff000000, white=0xffFFFFFF, // set more colors using Menu >  Tools > Color Selector
	   red=0xffFF0000, green=0xff00FF01,purple=0xffA000FF, blue=0xff0300FF, yellow=0xffFEFF00, cyan=0xff00FDFF, magenta=0xffFF00FB, grey=0xff5F5F5F, brown=0xffAF6407,
	   sand=0xffFCBA69, pink=0xffFF8EE7 ;
	// ************************************************************************ GRAPHICS 
	public void pen(int c, float w) {stroke(c); strokeWeight(w);}
	public void showDisk(float x, float y, float r) {ellipse(x,y,r*2,r*2);}
	
	// ************************************************************************ SAVING INDIVIDUAL IMAGES OF CANVAS 
	boolean snapPic=false;
	String PicturesOutputPath="data/PDFimages";
	int pictureCounter=0;
	public void snapPicture() {saveFrame("PICTURES/P"+nf(pictureCounter++,3)+".jpg"); }
	
	//************************ SAVING IMAGES for a MOVIE
	boolean filming=false;  // when true frames are captured in FRAMES for a movie
	int frameCounter=0;     // count of frames captured (used for naming the image files)
	boolean change=false;   // true when the user has presed a key or moved the mouse
	boolean animating=false; // must be set by application during animations to force frame capture
	/*
	To make a movie : 
	Press '~' to start filming, 
	act the movie or start an animation, 
	press '~' to pause/stop (you can restart to add frames)
	Then, from within your Processing sketch, 
	from the processing menu, select Tools > Movie Maker. 
	Click on Choose\u2026 Navigate to your Sketch Folder. 
	Select, but do not open, the FRAMES folder.
	Press Create Movie, 
	Select the parameters you want. 
	
	May not work for a large canvas!
	*/
	
	// ************************************************************************ TEXT 
	Boolean scribeText=true; // toggle for displaying of help text
	public void scribe(String S, float x, float y) {fill(0); text(S,x,y); noFill();} // writes on screen at (x,y) with current fill color
	public void scribeHeader(String S, int i) { text(S,10,20+i*20); noFill();} // writes black at line i
	public void scribeHeaderRight(String S) {fill(0); text(S,width-7.5f*S.length(),20); noFill();} // writes black on screen top, right-aligned
	public void scribeFooter(String S, int i) {fill(0); text(S,10,height-10-i*20); noFill();} // writes black on screen at line i from bottom
	public void scribeAtMouse(String S) {fill(0); text(S,mouseX,mouseY); noFill();} // writes on screen near mouse
	public void scribeMouseCoordinates() {fill(black); text("("+mouseX+","+mouseY+")",mouseX+7,mouseY+25); noFill();}
	public void displayHeader() { 
	    scribeHeader(title,0); 	   
	    }
	public void displayFooter() { // Displays help text at the bottom
	    scribeFooter(guide,1); 
	    scribeFooter(menu,0); 
	    }
	
	
	static public void main(String[] passedArgs) {
	    String[] appletArgs = new String[] { "SpiralFFTDesPKG.FFTSpiralDsn" };
	    if (passedArgs != null) {
	      PApplet.main(concat(appletArgs, passedArgs));
	    } else {
	      PApplet.main(appletArgs);
	    }
	  }
}
