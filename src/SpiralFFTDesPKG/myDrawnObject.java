package SpiralFFTDesPKG;
//Lots of stuff in here that is unnecessary for this project 
import java.util.ArrayList;
import java.util.Arrays;

import SpiralFFTDesPKG.FFTSpiralDsn.*;
import processing.core.PApplet;

public abstract class myDrawnObject {
	public static FFTSpiralDsn pa;

	public int[] clr;
	public float len;								//length of object
	protected static final int numReps = 4;				//default number of repetitions of subdivid/tuck/untuck
	
	protected pt[] origpts;							//originally drawn points making up this curve
	protected pt[] pts;								//points making up this curve
	protected float[] d_pts;						//array holding distance to each point from beginning
	
	//beautiful pts
	public cntlPt[] cntlPts;						//control points describing object, if used	
	
	protected float[] d_cntlPts;
	protected float cntl_len;
	
	public pt COV,									//center of verts
			COM;
	//boolean flags about object
	public boolean[] flags;							//init'ed to all false, then set specifically when appropriate
	//flag idxs
	public final int 
		//construction
			isClosed 		= 0,					//object is a closed poly
			isPosArea 		= 1,					//positive area: obj is oriented clockwise
			isMade 			= 2,					//whether or not the object is finished being drawn
			usesCntlPts 	= 3,					//if this curve is governed by control points (as opposed to just drawn freehand)
			firstDrawn		= 4,
	    //calculation
			doSweep 		= 5,					//whether or not the sweep should be -calculated- (used so that sweep is only made 1 time from draw for any config)
			reCalcPoints	= 6,					//recalculate points from cntl point radii - use if radius is changed on stroke from user input
			cntlWInvRad		= 7,					//whether the weight impact on cntl radius is inverse or direct - inverse has slow drawing be wider, direct has slow drawing be narrower
			useBallMorph	= 8,
			interpStroke	= 9,					//whether or not control-point based strokes are interpolating or not
		//display			 
			showWireFrame 	= 10,					//show wireframe version of sweep
			showSweep 		= 11,					//show the rotational solid construct	
			showCntlPnts 	= 12,					//show this object's cntl points
			vertNorms 		= 13,					//use vertex normals to shade curve
			drawNorms 		= 14,					//display normals for this object as small arrows
			drawCntlRad		= 15,
			useProcCurve	= 16;					//toggles whether we use straight lines in vertex building or processing's curve vertex			
	
	public final int numFlags = 17;					//always 1 more than last flag const
	
	//flags about type of operation that uses interpolation being done
	public int lnI_Typ,								//what interpolation type will this curve use for line operations (tuck, find pt at ratio of length, etc) 
				sbI_Typ,							//interp type for subdivision
				swI_Typ;							//what kind of interpolation will be used for this curve as it is swept around axis (if this is a closed sweep-poly)
	
	//flags about which interpolation type should be done
	public final int linear_int = 0,				//denotes linear interpolation 
					ball_int = 1;					//ball-morph interpolation, for sweep morph
	//add more here when we have more 
	
	public final int 	//point processing flags
		_subdivide		=0,
		_tuck			=1,
		_equaldist		=2,
		_resample		=3;

	//array of normals, tans for every control point,
	protected vec[] c_nAra, c_tAra;
	
	
	protected float[][] distFromIntAxis;			//keeps the distance to internal axis of stroke			
	protected static final int 
	//idxs in distFromAxis array
			_d = 0,									//index for dist from axis
			_t = 1,									//interpolant along rot axis for proj
			_a = 2,									//idx1 for pt in ara for rotational axis - proj lies between these two
			_b = 3;									//idx2 for pt in ara for rotational axis
	protected final int nAxisVals = 4;
	
	public myDrawnObject(FFTSpiralDsn _pa) {
		pa = _pa;
		initFlags();
		lnI_Typ = linear_int;
		swI_Typ = linear_int;
		sbI_Typ = linear_int;
		c_nAra = new vec[0];c_tAra = new vec[0];
		COM = pa.P();COV = pa.P();
		flags[firstDrawn] = true;
	}
	
	//initialize point referencing structs - using both for speed concerns.
	public void startDrawing(){	pts = new pt[0];len = 0; d_pts = new float[0]; flags[drawCntlRad] = false;}	
	public void addPt(pt p){
		ArrayList<pt> tmp = new ArrayList<pt>(Arrays.asList(pts));
		tmp.add(p);
		setPts(tmp);
	}//setPt
	
	//subdivide, tuck, respace, etc, cntlpts of this curve
	public void processCntlPts(int numPts, int numReps){
		float origLen = cntl_len;
		setCPts(procCntlPt(_subdivide, cntlPts, 2, origLen));										//makes 1 extra vert  equilspaced between each vert, to increase resolution of curve
		for(int i = 0; i < numReps; ++i){
			//setCPts(procCntlPt(_subdivide, cntlPts, 2, origLen));
			setCPts(procCntlPt(_tuck, cntlPts, .5f, origLen));
			setCPts(procCntlPt(_tuck, cntlPts, -.5f, origLen));
		}		//smooth curve - J4
		setCPts(procCntlPt(_equaldist, cntlPts, .5f, origLen));
		for(int i = 0; i < numReps; ++i){
			//setCPts(procCntlPt(_subdivide, cntlPts, 2, origLen));
			setCPts(procCntlPt(_tuck, cntlPts, .5f, origLen));
			setCPts(procCntlPt(_tuck, cntlPts, -.5f, origLen));
		}		//smooth curve - J4
		setCPts(procCntlPt(_resample, cntlPts, numPts, origLen));
	}			

	//subdivide, tuck, respace, resample, etc. pts of this curve
	public void processPts(pt[] _pts, int numPts, int numReps){
		setPts(procPts(_subdivide, _pts, 2, len, flags[isClosed]));										//makes 1 extra vert  equilspaced between each vert, to increase resolution of curve
		for(int i = 0; i < numReps; ++i){
			setPts(procPts(_subdivide, _pts, 2, len, flags[isClosed]));
			setPts(procPts(_tuck, _pts, .5f, len, flags[isClosed]));
			setPts(procPts(_tuck, _pts, -.5f, len, flags[isClosed]));
		}		//smooth curve - J4
		setPts(procPts(_equaldist, _pts, .5f, len, flags[isClosed]));
		for(int i = 0; i < numReps; ++i){
			setPts(procPts(_subdivide, _pts, 2, len, flags[isClosed]));
			setPts(procPts(_tuck, _pts, .5f, len, flags[isClosed]));
			setPts(procPts(_tuck, _pts, -.5f, len, flags[isClosed]));
		}		//smooth curve - J4
		setPts(procPts(_resample, _pts, numPts, len, flags[isClosed]));		
	}	
	
	
	//sets required info for points array - points and dist between pts, length, etc
	protected void setPts(ArrayList<pt> tmp){
		pts = tmp.toArray(new pt[0]);
		d_pts = getPtDist(pts, flags[isClosed]);	
		len=length(pts, flags[isClosed]);
	}//setPts	
	//make a new point interpolated between either 2 or 3 points in pts ara, described by # of idxs
	public pt makeNewPoint(pt[] _pts, int[] idxs, float s){	return _Interp(_pts[idxs[0]], s, (idxs.length == 2 ? _pts[idxs[1]] : _Interp(_pts[idxs[1]],.5f,_pts[idxs[2]], lnI_Typ)),lnI_Typ );	}
	
	/**
	 * process all points using passed algorithm on passed array of points - not all args are used by all algs.
	 * @param _typ type of point processing
	 * @param _pts array to be processed
	 * @param val quantity used by variou processing : subdivision-> # of new pts +1, tuck-> amt to be tucked,  resample-> # of new verts
	 * @param len length of segment described by points, including ends if closed
	 * @param wrap whether the point list wraps around or not
	 * @return arraylist of processed points
	 */
	public ArrayList<pt> procPts(int _typ, pt[] _pts, float val, float _len, boolean wrap){
		ArrayList<pt> tmp = new ArrayList<pt>(); // temporary array
		switch(_typ){
			case _subdivide	:{
			    for(int i = 0; i < _pts.length-1; ++i){tmp.add(_pts[i]); for(int j=1;j<val;++j){tmp.add(makeNewPoint(_pts,new int[]{i,i+1}, (j/(val))));}}
			    tmp.add(_pts[_pts.length-1]);				
			    return tmp;}
			case _tuck		:{
				if(wrap){tmp.add(makeNewPoint(_pts,new int[]{0,_pts.length-1,1}, val));} else {tmp.add(0,_pts[0]);}
			    for(int i = 1; i < _pts.length-1; ++i){	tmp.add(i,makeNewPoint(_pts,new int[]{i,i-1,i+1}, val));   }
		    	if(wrap){tmp.add(makeNewPoint(_pts,new int[]{_pts.length-1,_pts.length-2,0}, val));} else {tmp.add(_pts[_pts.length-1]);}			
		    	return tmp;}
			case _equaldist	:{
				float ratio = _len/(1.0f * _pts.length),curDist = 0;					 //new distance between each vertex, iterative dist travelled so far			 
				for(int i =0; i<_pts.length; ++i){tmp.add(at(curDist/_len));curDist+=ratio;}	
				tmp.add(_pts[_pts.length-1]);				
				return tmp;}	
			case _resample	:{
				float ratio = _pts.length/(1.0f * (val-1)),f;					//distance between each vertex		 
				int idx, newIdx=0;		
				for(float i = 0; i<_pts.length-1; i+=ratio){idx = (int)i;	f = i-idx;tmp.add(newIdx++,makeNewPoint(_pts,new int[]{idx,idx+1},f));}
				if(wrap) {
					if(pa.d(tmp.get(newIdx-1), tmp.get(0)) > ratio){	tmp.add(makeNewPoint(new pt[]{tmp.get(newIdx-1), tmp.get(0)},new int[]{0,1},.5f));}		//want to only add another point if last 2 points are further than ratio appart
				} else {		tmp.add(_pts[_pts.length-1]);}			//always add another point if open line/loop - want to preserve end point
				break;}	
			default :
		}
		
		return tmp;
	}
		
	//CONTROL POINT-RELATED FUNCTIONS
	//build essential orientation vectors for control points
	public void buildCntlFrameVecAras(){
		c_tAra = buildTangents(cntlPts, false);
		c_nAra = buildNormals(cntlPts,c_tAra);//rotate tans 90
	}//buildCntlFrameVecAras
		
	//sets required info for points array - points and dist between pts, length, etc
	protected void setCPts(ArrayList<cntlPt> tmp){
		cntlPts = tmp.toArray(new cntlPt[0]);
		d_cntlPts = getPtDist(cntlPts, false);	
		cntl_len=length(cntlPts, false);
	}//setPts	
	//make a new point interpolated between either 2 or 3 points in pts ara, described by # of idxs
	public cntlPt makeNewPoint(cntlPt[] _pts, int[] idxs, float s){	return _Interp(_pts[idxs[0]], s, (idxs.length == 2 ? _pts[idxs[1]] : _Interp(_pts[idxs[1]],.5f,_pts[idxs[2]], lnI_Typ)),lnI_Typ );	}
	/**
	 * process all points using passed algorithm on passed array of points - not all args are used by all algs.
	 * @param _typ type of point processing
	 * @param _pts array to be processed
	 * @param val quantity used by variou processing : subdivision-> # of new pts +1, tuck-> amt to be tucked,  resample-> # of new verts
	 * @param len length of segment described by points, including ends if closed
	 * @param wrap whether the point list wraps around or not
	 * @return arraylist of processed points
	 */	
	public ArrayList<cntlPt> procCntlPt(int _typ, cntlPt[] _pts, float val, float _len){
		ArrayList<cntlPt> tmp = new ArrayList<cntlPt>(); // temporary array
		switch(_typ){
			case _subdivide	:{
			    for(int i = 0; i < _pts.length-1; ++i){tmp.add(_pts[i]); for(int j=1;j<val;++j){tmp.add(makeNewPoint(_pts,new int[]{i,i+1}, (j/(val))));}}
			    tmp.add(_pts[_pts.length-1]);				
			    return tmp;}
			case _tuck		:{
				tmp.add(0,_pts[0]);//no wrap on control points, so no  need to check
			    for(int i = 1; i < _pts.length-1; ++i){	tmp.add(i,makeNewPoint(_pts,new int[]{i,i-1,i+1}, val));   }
		    	tmp.add(_pts[_pts.length-1]);			
		    	return tmp;}
			case _equaldist	:{
				float ratio = _len/(1.0f * _pts.length),curDist = 0;					 //new distance between each vertex, iterative dist travelled so far			 
				for(int i =0; i<_pts.length; ++i){tmp.add(at_C(curDist/_len, _pts));curDist+=ratio;}	
				tmp.add(_pts[_pts.length-1]);				
				return tmp;}	
			case _resample	:{
				float ratio = _pts.length/(1.0f * (val-1)),f;					//distance between each vertex		 
				int idx, newIdx=0;		
				for(float i = 0; i<_pts.length-1; i+=ratio){idx = (int)i;	f = i-idx;tmp.add(newIdx++,makeNewPoint(_pts,new int[]{idx,idx+1},f));}			
				tmp.add(_pts[_pts.length-1]);
				break;}	
			default :
		}
		return tmp;
	}	
	//end cntlpt related
	
	//normals, tangents, binormals at each point
	//build tangents first
	public vec[] buildTangents(pt[] _pts, boolean close){
		ArrayList<vec> tmp = new ArrayList<vec>();
		for(int i=0; i<_pts.length-1; ++i){tmp.add(pa.U(_pts[i], _pts[i+1]));}
		if(close){tmp.add(pa.U(_pts[_pts.length-1], _pts[0]));}
		else {tmp.add(pa.U(_pts[_pts.length-2], _pts[_pts.length-1]));}
		return tmp.toArray(new vec[0]);
	}	
	//must be called after tans are built -rotate tans 90
	public vec[] buildNormals(pt[] _pts, vec[] tAra){
		ArrayList<vec> tmp = new ArrayList<vec>();				
		for(int i =0; i<_pts.length; ++i){tmp.add(pa.R(tAra[i]));	}		//make normal the tangent rotated 90
		return tmp.toArray(new vec[0]);
	}

	//find location of center of verts
	public pt calcCOV(){pt C = pa.P();for(int i=0;i<pts.length;++i){C.add(pts[i]);}pt Ct = pa.P(1.0f/pts.length,C); COV=pa.P(Ct);return COV;}
	//find COV of passed verts
	public pt calcCOVOfAra(pt[] _pts){pt C = pa.P();for(int i=0;i<_pts.length;++i){C.add(_pts[i]);}pt Ct = pa.P(1.0f/_pts.length,C); return Ct;}

	
	/**
	 * return the interpolated vector between two pt's vectors given the adjacent idx's of two points in pts and the array of vectors
	 * @param idxA, idxB : 2 idxs in vec aras to be interped between
	 * @param s : interpolant
	 * @param vAra : array to find vecs in
	 * @return interpolated vector
	 */
	public vec getInterpVec(int idxA, float s, int idxB, vec[] vAra, int interpMech){return _Interp(vAra[idxA], s, vAra[idxB], interpMech);}
	public vec getUInterpVec(int idxA, float s, int idxB, vec[] vAra, int interpMech){return pa.U(_Interp(vAra[idxA], s, vAra[idxB], interpMech));}
	
	/**
	 * using arc length parameterisation this will return a point along the curve at a 
	 * particular fraction of the length of the curve (0,1 will return endpoints, .5 will return halfway along curve)
	 * @param t fraction of curve length we are interested in returning a point - should be 0-1
	 * @return point @t along curve
	 */
	public pt at(float t){return at(t,new float[1], len, pts, d_pts);}//put interpolant between adjacent axis points in s ara if needed
	public pt at(float t, float[] s){return at(t,s, len, pts, d_pts);}//put interpolant between adjacent axis points in s ara if needed
	public pt at(float t, float[] s, float _len, pt[] _pts, float[] _dpts){//call directly if wanting interpolant between adj axis points too
		if(t<0){PApplet.println("In at : t="+t+" needs to be [0,1]");return _pts[0];} else if (t>1){PApplet.println("In at : t="+t+" needs to be [0,1]");return _pts[_pts.length-1];}
		float dist = t * _len;
		for(int i=0; i<_dpts.length-1; ++i){										//built off d_pts so that it will get wrap for closed curve
			if(_dpts[i+1] >= dist){
				s[0] = ((dist-_dpts[i])/(_dpts[i+1]-_dpts[i]));					//needs to stay between 0 and 1 (since interpolation functions between pts will be 0-1 based), so normalize by distance d_pts[i]
				return makeNewPoint(_pts,new int[]{i,((i+1)%_pts.length)}, s[0]);		//put interpolant between adjacent axis points in s ara if needed		
			}					
		}		
		return _pts[0];
	}//at	
	
	public cntlPt at_C(float t, cntlPt[] _pts){float[] _dpts = this.getPtDist(_pts, false);float _len = this.length(_pts, false);return at_C(t,new float[1], _len, _pts, _dpts);}//put interpolant between adjacent axis points in s ara if needed
	public cntlPt at_C(float t, float[] s, float _len, cntlPt[] _pts, float[] _dpts){//call directly if wanting interpolant between adj axis points too
		if(t<0){PApplet.println("In at : t="+t+" needs to be [0,1]");return _pts[0];} else if (t>1){PApplet.println("In at : t="+t+" needs to be [0,1]");return _pts[_pts.length-1];}
		float dist = t * _len;
		for(int i=0; i<_dpts.length-1; ++i){										//built off d_pts so that it will get wrap for closed curve
			if(_dpts[i+1] >= dist){
				s[0] = ((dist-_dpts[i])/(_dpts[i+1]-_dpts[i]));					//needs to stay between 0 and 1 (since interpolation functions between pts will be 0-1 based), so normalize by distance d_pts[i]
				return makeNewPoint(_pts,new int[]{i,((i+1)%_pts.length)}, s[0]);		//put interpolant between adjacent axis points in s ara if needed		
			}			
		}		
		return new cntlPt(pa);
	}//at_C	
	

	
	/**
	 * this will conduct the appropriate interpolation on the two passed points, based on what type is set for the interpolation requested
	 * as more interpolation schemes are implemented, add cases
	 * @param A, B : points to be interpolated
	 * @param s : interpolant
	 * @param _typ : what is being interpolated - smoothing a line, sweeping the curve around the axis, etc.  for each _typ, a value should be set for this curve (i.e.smInterpTyp)
	 * @return : resultant point
	 */
	protected pt _Interp(pt A, float s, pt B, int _typ){
		switch (_typ){
			case linear_int : {	return pa.L(A, s, B);}
			//add more cases for different interpolation		
			default : {	return pa.L(A, s, B);}			//defaults to linear
		}	
	}//_Interp
	//same as above, but with vectors
	protected vec _Interp(vec A, float s, vec B, int _typ){
		switch (_typ){
			case linear_int : {	return pa.L(A, s, B);}
			//add more cases for different interpolation		
			default : {	return pa.L(A, s, B);}			//defaults to linear
		}	
	}//_Interp
	//same as above but with doubles
	protected float _Interp(float A, float s, float B, int _typ){
		switch (_typ){
			case linear_int : {	return (1-s)*A + (s*B);}
			//add more cases for different interpolation		
			default : {	return (1-s)*A + (s*B);}			//defaults to linear
		}	
	}//_Interp
	protected cntlPt _Interp(cntlPt A, float s, cntlPt B, int _typ){
		switch (_typ){
			case linear_int : {	return cntlPt.L(A, s, B);}
			//add more cases for different interpolation		
			default : {	return cntlPt.L(A, s, B);}			//defaults to linear
		}	
	}//_Interp	

	//draw currently selected point
	public void drawSelPoint(int i ){
		pa.pushMatrix();		pa.pushStyle();
		pa.stroke(255,255,0,255);
		if(flags[usesCntlPts]){cntlPts[i].show(3);} else {pts[i].show(3);}
		pa.popStyle();		pa.popMatrix();
	}
	

	public abstract void rebuildPolyPts();

	//makes a copy of the points in order
	public pt[] cpyPoints(pt[] _pts){pt[] tmp = new pt[_pts.length]; for(int i=0; i<_pts.length; ++i){	tmp[i]=pa.P(_pts[i]);}	return tmp;}//cpyPoints
	//makes a copy of the points in order
	public cntlPt[] cpyPoints(cntlPt[] _pts){cntlPt[] tmp = new cntlPt[_pts.length]; for(int i=0; i<_pts.length; ++i){	tmp[i]=new cntlPt(_pts[i]);}	return tmp;}//cpyPoints

	//move points by the passed vector 
	public pt[] movePoints(vec move, pt[] _pts){for(int i =0; i<_pts.length; ++i){	_pts[i].add(move);	}	return _pts;}	
	//move points by the passed vector 
	public cntlPt[] movePoints(vec move, cntlPt[] _pts){for(int i =0; i<_pts.length; ++i){	_pts[i].add(move);	}	return _pts;}	
	//flip the passed points and move them based on the displacement from the passed movement vector
//	public pt[] flipPtsAndMove(myDrawnObject _obj, pt[] _pts, vec move,  vec covAxis){
//		pt[] tmp = movePoints(move, cpyPoints(_pts));
//		return tmp;
//	}//flipPtsAndMove
//	//flip the passed cntlpts and move them based on the displacement from the passed movement vector
//	public cntlPt[] flipPtsAndMove(myDrawnObject _obj, cntlPt[] _pts, vec move,  vec covAxis, boolean reverse){
//		cntlPt[] tmp = movePoints(move, (reverse ? rotPtsAroundCOV(_pts, PApplet.PI, _obj.COV, pa.U(_obj.canvasNorm), covAxis) : cpyPoints(_pts)));
//		return tmp;
//	}//flipPtsAndMove
	//set this object's points to be passed points, for copying
	public void setPtsToArrayList(cntlPt[] _pts){ArrayList<cntlPt> tmp = new ArrayList<cntlPt>(Arrays.asList(_pts));setCPts(tmp);}	
	//set this object's points to be passed points, for copying
	public void setPtsToArrayList(pt[] _pts){ArrayList<pt> tmp = new ArrayList<pt>(Arrays.asList(_pts));setPts(tmp);}	

	
	//rotate points around axis that is xprod of canvas norm and the lstroke cov to the rstroke cov vec at stroke cov.
	//should make mirror image of pts
	public pt[] rotPtsAroundCOV(pt[] _pts, float angle, pt cov, vec _canvasNorm, vec _covNorm){//need to pass canvas norm since it cannot be static
		ArrayList<pt> tmp = new ArrayList<pt>();//				res.get(sl)[i] = pa.R(old, sliceA, canvasNorm, bv, ptOnAxis);
		for(int i=0; i<_pts.length; ++i){tmp.add(pa.R(_pts[i], angle, pa.V(_canvasNorm), _covNorm, cov));}
		return tmp.toArray(new pt[0]);			
	}//	
	//rotate points around axis that is xprod of canvas norm and the lstroke cov to the rstroke cov vec at stroke cov.
	//should make mirror image of pts
	public cntlPt[] rotPtsAroundCOV(cntlPt[] _pts, float angle, pt cov, vec _canvasNorm, vec _covNorm){//need to pass canvas norm since it cannot be static
		ArrayList<cntlPt> tmp = new ArrayList<cntlPt>();//				res.get(sl)[i] = pa.R(old, sliceA, canvasNorm, bv, ptOnAxis);
		for(int i=0; i<_pts.length; ++i){tmp.add(pa.R(_pts[i], angle, pa.V(_canvasNorm), _covNorm, cov));}
		return tmp.toArray(new cntlPt[0]);			
	}//	
	//finds index of point with largest projection on passed vector in passed pt ara
	protected int findLargestProjection(vec v, pt c, pt[] _pts){
		float prjLen = -1, d;
		int res = -1;
		for(int i=0; i<_pts.length; ++i){d = pa.dot(v,pa.V(c, _pts[i]));if(d > prjLen){prjLen = d;res = i;}}	
		return res;
	}//findLargestProjection : largest projection on passed vector in passed pt ara		
	
	public abstract int findClosestPt(pt p, float[] d);
	//finds closest point to p in sPts - put dist in d
	protected final int findClosestPt(pt p, float[] d, pt[] _pts){
		int res = -1;
		float mindist = 99999999, _d;
		for(int i=0; i<_pts.length; ++i){_d = pa.d(p,_pts[i]);if(_d < mindist){mindist = _d; d[0]=_d;res = i;}}	
		return res;
	}
	//reorder verts so that they start at newStIdx - only done for closed  meshes, so should wrap around
	//made static so doesn't have to be used with pts array - should be called by implementing class
	protected ArrayList<pt> reorderPts(int newStIdx, pt[] _pts){
		ArrayList<pt> tmp = new ArrayList<pt>(Arrays.asList(_pts));
		for(int i = 0; i<_pts.length; ++i){tmp.set(i, _pts[(i+newStIdx)%_pts.length]);	}
		return tmp;
	}		
	
	public abstract void remakeDrawnTraj(boolean useVels);

	//run at the end of drawing a curve - will set appropriate flags, execute smoothing, subdivision, resampling, etc
	public abstract void finalizeDrawing(boolean procPts);	

	public void drawMe() {
		pa.pushMatrix();
		pa.pushStyle();
			pa.fill(clr[0],clr[1],clr[2],50);
			pa.stroke(0,0,0,255);
			pa.strokeWeight(1);
//			if(flags[useProcCurve]){pa.show(pts);} 
//			else {			
				pa.curve(pts);
				//}
		pa.popStyle();			
		pa.popMatrix();
//		if(flags[drawNorms] && (nAra != null)&& (tAra != null)&& (bAra != null)){drawNorms(pts, nAra,tAra,bAra,20);}
//		drawCOV();
	}
	public abstract void dragPicked(vec disp, int idx);
	public abstract void dragAll(vec disp);			//move COV to location pointed at by mouse, and move all points by same displacement
	//drag point represented by passed idx in passed array - either point or cntl point
	protected void dragPicked(vec dispVec, int idx, pt[] _pts) {if((-1 != idx) && (_pts[idx] != null)){_pts[idx].add(dispVec);flags[reCalcPoints]=true;}}
	protected void dragAll(vec dispVec, pt[] _pts){if((_pts == null)||(_pts.length == 0)){return;}for(int i=0;i<_pts.length;++i){_pts[i].add(dispVec);}flags[reCalcPoints]=true;}
	
	
	//returns array of distances to each point from beginning - needs to retain dist from last vert to first if closed
	//public final float[] getPtDist(){float[] res = new float[pts.length+1];res[0]=0;for(int i=1; i<pts.length; ++i){res[i] = res[i-1] + pa.d(pts[i-1],pts[i]);}if(flags[isClosed]){res[pts.length] = res[pts.length-1] +  pa.d(pts[pts.length-1],pts[0]);}return res;}
	//public final float[] getPtDist(){return getPtDist(pts, flags[isClosed]);}
	//returns array of distances to each point from beginning - needs to retain dist from last vert to first if closed
	public final float[] getPtDist(pt[] _pts, boolean wrap){
		float[] res = new float[_pts.length+1];
		res[0]=0;
		for(int i=1; i<_pts.length; ++i){
			//System.out.println("i : "+i);
			res[i] = res[i-1] + pa.d(_pts[i-1],_pts[i]);
			}
		if(wrap){
			//System.out.println("wrap");
			res[_pts.length] = res[_pts.length-1] +  pa.d(_pts[_pts.length-1],_pts[0]);
		} else {
			//System.out.println("no wrap");
			
			res[_pts.length] = 0;
		}
		
		return res;}
	//returns length of curve, including endpoint if closed
	//public final float length(){return length(pts, flags[isClosed]);}//{float res = 0;for(int i =0; i<pts.length-1; ++i){res += pa.d(pts[i],pts[i+1]);}if(flags[isClosed]){res+=pa.d(pts[pts.length-1],pts[0]);}return res;}
	public final float length(pt[] _pts, boolean closed){float res = 0;for(int i =0; i<_pts.length-1; ++i){res += pa.d(_pts[i],_pts[i+1]);}if(closed){res+=pa.d(_pts[_pts.length-1],_pts[0]);}return res;}


	public void drawCOV(){		if(COV == null) {return;}		pa.pushMatrix();		pa.pushStyle();	pa.stroke(255,0,255,255);		COV.show(3);		pa.popStyle();		pa.popMatrix();	}
	//drawCntlRad
	public pt getPt(int i){return pts[i];}
	
	public void initFlags(){flags = new boolean[numFlags]; for(int i=0;i<numFlags;++i){flags[i]=false;}}
	public String toString(){
		String res = "#pts : "+pts.length+" len : "+ len ;
		return res;
	}
	
	public void toggleHatch() {
	    
	}
}

