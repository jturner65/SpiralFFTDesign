package SpiralFFTDesPKG;
import processing.core.PApplet;
import processing.core.PConstants;

import java.util.ArrayList;
import java.util.Arrays;

import SpiralFFTDesPKG.FFTSpiralDsn.*;


//All user interaction with this curve is through the cntlpoints.
public class myVariStroke extends myDrawnObject {

	protected final int numVerts = 200;							

	public int offsetType;						//1:Q-bspline w/normal offset, 2:Q-bspline w/ball offset, 3:Q-bspline w/radial offset
	public myOffset _offset;					//offset used by this stroke to calculate poly loop
	public final int numIntCntlPts = 200, numCntlPts = 6;			//# of control points to use to draw line
	private boolean ptsDerived;					//whether or not the points making up the loop of this stroke have been derived yet

	
	public pt[] drawnCntlPts;						//control point locations - start point +  vel vector (scaled tangent), repeat
	public float[] d_drwnCntlPts;					//distance between each drawnCntlPts points
	public float drwnCntlLen;						//len of arc of drawnCntlPts pts

	public float[] cntlPtIntrps;				//interpolants for each control point as part of total, based upon radius of controlpoint - larger will be bigger. 
												//for each cntl pt, this value will be cntl point rad / total control point rads.
	//interpolated drawn curve, weighted by drawn speed
	public pt[] interpCntlPts;					//interpolated control point locations - start point +  vel vector (scaled tangent), repeat
	public float[] d_interpPts;					//distance between interpolated points
	public float interpLen;						//len of interp pts
	
	public myVariStroke(FFTSpiralDsn _pa) {
		super(_pa);
		flags[isClosed] = false;		
		clr =new int[]{255,0,0,255}; //red
		flags[drawCntlRad] = true;
	    cntlPts = new cntlPt[0];
	    interpCntlPts = new pt[0];
		_offset = new myNormOffset(_pa);
		flags[usesCntlPts] = true;
		flags[interpStroke] = true;
		ptsDerived = false;
		flags[cntlWInvRad] = false;			//whether slow drawing makes rad larger or smaller
	}

	//as drawing, add points to -cntlPoints-, not pts array.
	public void addPt(pt p){
		ArrayList<cntlPt> tmp = new ArrayList<cntlPt>(Arrays.asList(cntlPts));
		int i = tmp.size()-1;
		if(i > 0 ){tmp.get(i).w = calcCntlWeight(p,tmp.get(i),tmp.get(i-1));}//previous point's weight 
		cntlPt tmpPt = new cntlPt(pa, p);
		tmp.add(tmpPt);
		setCPts(tmp);
	}//
	
	public void finalizeDrawing(boolean procPts){
		//by here we have the drawn points from user input
		//we want use the offset to determine the actual points of the curve
		//we want to put this in a function so that any changes to the 
		//cntlpoints can be cascaded down to the actual loop		
		buildPointsUsingOffset(procPts, numReps);
		//calculate line points from control points
		//find loop around stroke line by cntl points' radii
		//once loop is built, treat as poly loop
		flags[isMade] = true;
		flags[drawCntlRad] = false;
		//build array of weights from 
		buildInterpAra();
	}//finalize
	
	//build pts array using cntlpoints and offset chosen
	public void buildPointsUsingOffset(boolean procPts, int repCnt){
		if(procPts){
		    finalizeCntlW();
		    for(int i=0;i<cntlPts.length;++i){cntlPts[i].calcRadFromWeight(cntl_len/cntlPts.length, flags[cntlWInvRad]);}           //sets all radii based on weights
		    processCntlPts(flags[interpStroke] ? numIntCntlPts : numCntlPts, repCnt);
	    }
		buildCntlFrameVecAras();
		buildPtsFromCntlPts();
	}
	//build interpolants based upon weight of each cntl point.  these can be used to determine how far the edge moves (velocity) and how much it rotates per frame
	//when this is done, we have multipliers to apply to each tangent vector to determine displacement for each of the frames
	public void buildInterpAra(){
		cntlPtIntrps = new float[numIntCntlPts];			//interpolant from each control point - weights
		
		float[] cptInterps = new float[numIntCntlPts];
		interpCntlPts = new pt[numIntCntlPts];
		drawnCntlPts = new pt[numIntCntlPts];
		float sumWts = 0;
		for(int i=0;i<cntlPts.length;++i){sumWts += cntlPts[i].w;cptInterps[i] = cntlPts[i].w;cntlPtIntrps[i]=sumWts;}
		//for(int i=0;i<cntlPts.length;++i){sumWts += cntlPts[i].w;cntlPtIntrps[i]=cntlPts[i].w;}
		//System.out.println("total weight = " + sumWts);
		for(int i=0;i<cntlPts.length;++i){cntlPtIntrps[i]/=sumWts;cptInterps[i]/=sumWts;}
		//smooth interpolants now
		cntlPtIntrps = dualFloats(cntlPtIntrps);
		cptInterps = dualFloats(cptInterps);
		
		interpCntlPts[0] = pa.P(cntlPts[0]);			//set first point
		drawnCntlPts[0] = pa.P(cntlPts[0]);
		float distStToEnd = pa.d(cntlPts[0], cntlPts[cntlPts.length-1]);			//distance from first to last point
		
		for(int i=1;i<cntlPts.length;++i){
			interpCntlPts[i] = pa.P(interpCntlPts[i-1],pa.W(distStToEnd * cptInterps[i], c_tAra[i]));		
			drawnCntlPts[i] = pa.P(cntlPts[i]);
		}	
		for(int i= cntlPts.length; i<numIntCntlPts; ++i){
			interpCntlPts[i] = pa.P(interpCntlPts[i-1],pa.W(distStToEnd * cptInterps[i], c_tAra[c_tAra.length-1]));		
			drawnCntlPts[i] = pa.P(cntlPts[cntlPts.length-1]);
			
		}
		//System.out.println("# interpCntlPts : "+ interpCntlPts.length + " cntlPts.length : " + cntlPts.length);		
		
		d_interpPts = getPtDist(interpCntlPts, false);	
		interpLen = length(interpCntlPts, false);
		d_drwnCntlPts = getPtDist(drawnCntlPts, false);	
		drwnCntlLen = length(drawnCntlPts, false);		
		//smooth/equi-space interpolated cntl points
		processInterpPts( numIntCntlPts, numReps);
		//for(int i=0;i<cntlPts.length;++i){System.out.println("Tot weight ["+i+"] = "+cntlPtIntrps[i]+ "\tweight at ["+i+"] = " + cptInterps[i] +  "\tCNTLPT : "+cntlPts[i].toString()+"|\tInterp Point : "+interpCntlPts[i].toString());}
	}	
	//subdivide, tuck, respace, resample, etc. pts of this curve
	public void processInterpPts(int numPts, int numReps){
		setInterpPts(procInterpPts(_subdivide, interpCntlPts, 2, interpLen));										//makes 1 extra vert  equilspaced between each vert, to increase resolution of curve
		for(int i = 0; i < numReps; ++i){
			//setInterpPts(procInterpPts(_subdivide, interpCntlPts, 2, interpLen));
			setInterpPts(procInterpPts(_tuck, interpCntlPts, .5f, interpLen));
			setInterpPts(procInterpPts(_tuck, interpCntlPts, -.5f, interpLen));
		}		//smooth curve - J4
		setInterpPts(procInterpPts(_equaldist, interpCntlPts, .5f, interpLen));
		for(int i = 0; i < numReps; ++i){
			//setInterpPts(procInterpPts(_subdivide, interpCntlPts, 2, interpLen));
			setInterpPts(procInterpPts(_tuck, interpCntlPts, .5f, interpLen));
			setInterpPts(procInterpPts(_tuck, interpCntlPts, -.5f, interpLen));
		}		//smooth curve - J4
		setInterpPts(procInterpPts(_resample, interpCntlPts, numPts, interpLen));	
	}	
	
	//return appropriate ara of points based on using velocities or not
	public pt[] getDrawnPtAra(boolean useVels){	return (useVels ? interpCntlPts : cntlPts);}
	
	//remake drawn trajectory after edit
	public void remakeDrawnTraj(boolean useVels){
		if(useVels){
			for(int i = 0; i < 10; ++i){
				//setInterpPts(procInterpPts(_subdivide, interpCntlPts, 2, interpLen));
				setInterpPts(procInterpPts(_tuck, interpCntlPts, .5f, interpLen));
				setInterpPts(procInterpPts(_tuck, interpCntlPts, -.49f, interpLen));
			}		//smooth curve - J4
			setInterpPts(procInterpPts(_equaldist, interpCntlPts, .5f, interpLen));
			setInterpPts(procInterpPts(_resample, interpCntlPts, numIntCntlPts, interpLen));	
		} else {
			for(int i = 0; i < 10; ++i){
				//setInterpPts(procInterpPts(_subdivide, interpCntlPts, 2, interpLen));
				setCPts(procCntlPt(_tuck, cntlPts, .5f, cntlPts.length));
				setCPts(procCntlPt(_tuck, cntlPts, -.49f, cntlPts.length));
			}		//smooth curve - J4
			setCPts(procCntlPt(_equaldist, cntlPts, .5f, cntlPts.length));
			setCPts(procCntlPt(_resample, cntlPts, numIntCntlPts, cntlPts.length));
		}
		System.out.println("# interpCntlPts : "+ interpCntlPts.length + " cntlPts.length : " + cntlPts.length);	
	}//remakeDrawnTraj
	
	
	public void handleMouseDrag(pt mse, vec dispVec, int drawnTrajPickedIdx){		
		pt[] pts = getDrawnPtAra(pa.useDrawnVels);
		int minBnd = PApplet.max(drawnTrajPickedIdx - pa.drawnTrajEditWidth, 0),
			maxBnd = PApplet.min(drawnTrajPickedIdx + pa.drawnTrajEditWidth, pts.length-1);		
		//System.out.println("drag in drag zone inside disp calc -> idx bounds : " + minBnd + " | " + maxBnd);
		float modAmt, invdistLow = 1.0f/(drawnTrajPickedIdx - minBnd), invdistHigh = 1.0f/(maxBnd - pa.drawnTrajPickedIdx);
		for(int idx = minBnd; idx < maxBnd; ++idx){
			float divMultVal = (idx > pa.drawnTrajPickedIdx) ? invdistHigh:invdistLow;
			modAmt = pa.trajDragScaleAmt* PApplet.cos((idx-drawnTrajPickedIdx) * PConstants.HALF_PI * divMultVal);//trajDragScaleAmt/abs(1 + (idx-drawnTrajPickedIdx));
			//modAmt *= modAmt;
			pts[idx].add(pa.W(modAmt,dispVec));
		}
		
	}
	//sets required info for points array - points and dist between pts, length, etc
	protected void setInterpPts(ArrayList<pt> tmp){
		interpCntlPts = tmp.toArray(new pt[0]);
		d_interpPts = getPtDist(interpCntlPts, false);	
		interpLen=length(interpCntlPts, false);
	}//setPts	
	//tuck untuck float values
	public float[] dualFloats(float[] src){
		float[] res = new float[src.length],res1 = new float[src.length];
		res1[0]=src[0];
		res1[src.length-1]=src[src.length-1];
		for(int i=1; i<src.length-1;++i){
			res1[i]=_Interp(src[i],.5f,_Interp(src[i-1],.5f,src[i+1],lnI_Typ),lnI_Typ);
		}
		res[0]=res1[0];
		res[src.length-1]=res1[src.length-1];
		for(int i=1; i<res1.length-1;++i){
			res[i]=_Interp(res1[i],-.5f,_Interp(res1[i-1],.5f,res1[i+1],lnI_Typ),lnI_Typ);
		}			
		return res;		
	}
	
	/**
	 * process all points using passed algorithm on passed array of points - not all args are used by all algs.
	 * @param _typ type of point processing
	 * @param _pts array to be processed
	 * @param val quantity used by variou processing : subdivision-> # of new pts +1, tuck-> amt to be tucked,  resample-> # of new verts
	 * @param len length of segment described by points, including ends if closed
	 * @return arraylist of processed points
	 */
	public ArrayList<pt> procInterpPts(int _typ, pt[] _pts, float val, float _len){
		ArrayList<pt> tmp = new ArrayList<pt>(); // temporary array
		switch(_typ){
			case _subdivide	:{
			    for(int i = 0; i < _pts.length-1; ++i){tmp.add(_pts[i]); for(int j=1;j<val;++j){tmp.add(makeNewPoint(_pts,new int[]{i,i+1}, (j/(val))));}}
			    tmp.add(_pts[_pts.length-1]);				
			    return tmp;}
			case _tuck		:{
				tmp.add(0,_pts[0]);
			    for(int i = 1; i < _pts.length-1; ++i){	tmp.add(i,makeNewPoint(_pts,new int[]{i,i-1,i+1}, val));   }
		    	tmp.add(_pts[_pts.length-1]);		
		    	return tmp;}
			case _equaldist	:{
				float ratio = _len/(1.0f * _pts.length),curDist = 0;					 //new distance between each vertex, iterative dist travelled so far			 
				for(int i =0; i<_pts.length; ++i){tmp.add(at_I(curDist/_len));curDist+=ratio;}	
				tmp.add(_pts[_pts.length-1]);				
				return tmp;}	
			case _resample	:{
				float ratio = _pts.length/(1.0f * (val-1)),f;					//distance between each vertex		 
				int idx, newIdx=0;		
				for(float i = 0; i<_pts.length-1; i+=ratio){idx = (int)i;	f = i-idx;tmp.add(newIdx++,makeNewPoint(_pts,new int[]{idx,idx+1},f));}
				tmp.add(_pts[_pts.length-1]);			//always add another point if open line/loop - want to preserve end point
				break;}	
			default :
		}
		
		return tmp;
	}

	public pt at_I(float t){return at(t,new float[1], interpLen, interpCntlPts, d_interpPts);}//put interpolant between adjacent points in s ara if needed
	public pt at_I(float t, float[] s){	return at(t,s, interpLen, interpCntlPts, d_interpPts);}//put interpolant between adjacent points in s ara if needed
	
	
	private void buildPtsFromCntlPts(){
		ArrayList<pt> tmp =  _offset.calcOffset(cntlPts, c_nAra, c_tAra) ;
		pts = tmp.toArray(new pt[0]);		
		ptsDerived = true;		
	}//buildPtsFromCntlPts

	//calculate initial weights for last few points of drawn cntlpoint line
	private void finalizeCntlW(){
		if(cntlPts.length == 0){return ;}
		cntlPts[0].w = cntlPts.length < 2 ? 1 : cntlPts[1].w;
		if(cntlPts.length == 1){return ;}
		cntlPts[cntlPts.length-2].w = cntlPts.length < 3 ? cntlPts[0].w : calcCntlWeight(cntlPts[cntlPts.length-3],cntlPts[cntlPts.length-2],cntlPts[cntlPts.length-3]);
		cntlPts[cntlPts.length-1].w = cntlPts.length < 2 ? cntlPts[0].w : cntlPts[cntlPts.length-2].w;
	}
	
	//build poly loop points using offsets from control points radii
	public void rebuildPolyPts(){
		flags[reCalcPoints]=false;
		//buildPointsUsingOffset(true,1);
		buildPtsFromCntlPts();
		flags[isMade] = true;
	}
	
	public void drawMe(boolean useDrawnVels){
		pa.pushMatrix();
		pa.pushStyle();
			pa.fill(clr[0],clr[1],clr[2],clr[3]);
			pa.stroke(0,0,0,255);
			pa.strokeWeight(1);
        	if(useDrawnVels){
        		int clrInt = 0;
    			for(int i = 0; i < interpCntlPts.length; ++i){
    	        	clrInt = (int)(i/(1.0f * interpCntlPts.length) * 255.0f);
    	            pa.fill(clrInt,255,(255 - clrInt),255);  
    	            pa.stroke(clrInt,255,(255 - clrInt),255); 
    				this.interpCntlPts[i].show(2);
    				if(flags[drawCntlRad]){pa.circle(this.interpCntlPts[i], this.cntlPts[i].r,this.c_nAra[i], this.c_tAra[i],20);}
    			}
        	} else {			
				for(int i = 0; i < cntlPts.length; ++i){
					cntlPts[i].drawMe(ptsDerived ? 0 : 1);
				}
				if(flags[drawCntlRad]){this._offset.drawCntlPts(this.cntlPts, this.c_nAra, this.c_tAra, ptsDerived);}
        	}
		pa.popStyle();			
		pa.popMatrix();		
	}//
	
	public pt[] moveVelCurveToEndPoints(pt startPt,pt endPt){
		pt[] destCurve = new pt[interpCntlPts.length];

		if(interpCntlPts.length == 0){return destCurve;}
		//drawn curve params
		pt origin = interpCntlPts[0];
		pt end = interpCntlPts[interpCntlPts.length - 1];
		
		//edge params		
		vec drawnAxis = pa.V(origin, end);
		vec edgeAxis =  pa.V(startPt, endPt);		//angle between these two is the angle to rotate everyone
		
		//transformation params
		vec dispToStart = pa.V(origin, startPt);			//displacement vector between start of drawn curve and edge 1.
		float alpha = -pa.angle(drawnAxis,edgeAxis);			//angle to rotate everyone
		float scaleRatio = pa.n(edgeAxis)/pa.n(drawnAxis);	//ratio of distance from start to finish of drawn traj to distance between edges - multiply all elements in drawn traj by this
		
		//displace to align with start
		destCurve = movePoints(dispToStart, interpCntlPts);

		//displace every point to be scaled distance from start of curve equivalent to scale of edge distances to drawn curve
		for(int ptItr = 1; ptItr < interpCntlPts.length ; ++ptItr){
			destCurve[ptItr].set(pa.P(destCurve[0],scaleRatio, pa.V(destCurve[0],destCurve[ptItr])));//start point displaced by scaleRatio * vector from start to original location of pt
		}
		//rotate every point around destCurve[0] by alpha
		for(int ptItr = 1; ptItr < interpCntlPts.length ; ++ptItr){
			destCurve[ptItr].rotate(alpha,destCurve[0]);
		}
		interpCntlPts = destCurve;
		return destCurve;
	}//fitCurveBetweenJT
	
	public cntlPt[] moveCntlCurveToEndPoints(pt startPt,pt endPt){
		cntlPt[] destCurve = new cntlPt[cntlPts.length];
		
		//drawn curve params
		if(cntlPts.length == 0){return destCurve;}
		pt origin = cntlPts[0];
		pt end = cntlPts[cntlPts.length - 1];
		
		//edge params		
		vec drawnAxis = pa.V(origin, end);
		vec edgeAxis =  pa.V(startPt, endPt);		//angle between these two is the angle to rotate everyone
		
		//transformation params
		vec dispToStart = pa.V(origin, startPt);			//displacement vector between start of drawn curve and edge 1.
		float alpha = -pa.angle(drawnAxis,edgeAxis);			//angle to rotate everyone
		float scaleRatio = pa.n(edgeAxis)/pa.n(drawnAxis);	//ratio of distance from start to finish of drawn traj to distance between edges - multiply all elements in drawn traj by this
		
		//displace to align with start
		destCurve = movePoints(dispToStart, cntlPts);

		//displace every point to be scaled distance from start of curve equivalent to scale of edge distances to drawn curve
		for(int ptItr = 1; ptItr < cntlPts.length ; ++ptItr){
			destCurve[ptItr].set(pa.P(destCurve[0],scaleRatio, pa.V(destCurve[0],destCurve[ptItr])));//start point displaced by scaleRatio * vector from start to original location of pt
		}
		//rotate every point around destCurve[0] by alpha
		for(int ptItr = 1; ptItr < cntlPts.length ; ++ptItr){
			destCurve[ptItr].rotate(alpha,destCurve[0]);
		}
		cntlPts = destCurve;
		return destCurve;
	}//fitCurveBetweenJT
	
	//move points by the passed vector 
	public cntlPt[] movePoints(vec move, cntlPt[] _pts){for(int i =0; i<_pts.length; ++i){	_pts[i].add(move);	}	return _pts;}
	
	//calculate the weight of each point by determining the distance from its two neighbors - radius is inversely proportional to weight
	public float calcCntlWeight(pt a, pt p, pt b){	return (pa.d(a,p) + pa.d(p,b));}
	
	//override this for cntrl-point driven constructs
	public int findClosestPt(pt p, float[] d){	
		return findClosestPt(p, d,cntlPts);
	}

	//drag point represented by passed idx in passed array - either point or cntl point
	public void dragPicked(vec disp, int idx) {dragPicked(disp, idx, cntlPts);}
	//drag all points by finding displacement to mouse for COV and moving all points by same dislacement
	public void dragAll(vec disp) {dragAll(disp, cntlPts);}	

	
	public String toString(){
		String res = "Interpolating Spline Stroke : Offset calc : "+ _offset;
		res += super.toString();
		return res;	
	}
	
}//myVariStroke
