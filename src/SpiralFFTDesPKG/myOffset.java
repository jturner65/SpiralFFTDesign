package SpiralFFTDesPKG;
import java.util.*;

import SpiralFFTDesPKG.FFTSpiralDsn.*;
import processing.core.PConstants;

//class to hold functionality to calculate offset "sidewalks"

public abstract class myOffset {
	public static FFTSpiralDsn pa;
	public int ID;
	public static int IDcnt = 0;
	public String name;
	public int capSize = 20;
	public boolean endCaps;

	public myOffset(FFTSpiralDsn _pa, boolean _ec){
		pa = _pa;
		ID = IDcnt++;
		endCaps = _ec;
	}	
	
	public myOffset(FFTSpiralDsn _pa){this(_pa, true);}
	
	/**
	 * calculate the offset points for the drawn stroke line contained in _obj
	 * @param _obj drawn stroke to build offset pts from
	 */
	public abstract ArrayList<pt> calcOffset(cntlPt[] cntlPts, vec[] nAra, vec[] tAra);				
	public abstract void drawCntlPts(cntlPt[] pts, vec[] nAra, vec[] tAra, boolean derived);
	
	/**
	 * build an array of points that sweeps around c clockwise in plane of norm and tan, with starting radius c.r * norm
	 * @param c control point
	 * @param norm normal of point (binormal in world frame, this is direction of offset)
	 * @param tan tangent at point
	 * @return pt array of sequence of points in an arc for an endcap
	 */
	public ArrayList<pt> buildCapPts (cntlPt c, vec norm, vec tan, float mult){
		ArrayList<pt> tmp = new ArrayList<pt>();
		float angle = PConstants.PI/(1.0f*capSize), sliceA = angle;			//10 slices
		tmp.add(pa.P((pt)c, mult * -c.r, norm));
		for(int i=1;i<capSize-1;++i){	tmp.add(pa.R(tmp.get(i-1), sliceA, norm, tan, c));}
		tmp.add(pa.P((pt)c, mult * c.r, norm));
		return tmp;	
	}//buildCapPts
	
	public String toString(){
		String res = name + "Offset ID : "+ID;		
		return res;
	}
}//myOffset

/**
 * calculates normal offset - distance r, normal from stroke line
 * @author john
 *
 */
//make other classes to use different offset mechanism
class myNormOffset extends myOffset{
	myNormOffset(FFTSpiralDsn _pa){super(_pa); name = "Normal offset";}

	@Override
	public  ArrayList<pt> calcOffset(cntlPt[] cntlPts, vec[] nAra, vec[] tAra) {
		if(nAra.length != cntlPts.length){return  new ArrayList<pt>();}	
		ArrayList<pt> tmp = new ArrayList<pt>();
		int numCptsM1 = cntlPts.length-1;		
		//start at first point and build endcap
		if(endCaps){tmp.addAll(buildCapPts(cntlPts[0], nAra[0], tAra[0], 1));}
		for(int i = 0; i<cntlPts.length;++i){	tmp.add(pa.P((pt)cntlPts[i], cntlPts[i].r, nAra[i]));}//add cntl point + rad offset from norm
		//build endcap on last cntlpoint
		if(endCaps){tmp.addAll(buildCapPts(cntlPts[numCptsM1], nAra[numCptsM1], tAra[numCptsM1],-1));}
		for(int i = numCptsM1; i>=0;--i){	tmp.add(pa.P((pt)cntlPts[i], -cntlPts[i].r, nAra[i]));}//add cntl point + rad offset from norm negated, in backwards order, so all points are added properly
		return tmp;
	}
	
	public String toString(){
		String res = name +super.toString();		
		return res;
	}
	
    @Override
//    public void drawCntlPts(cntlPt[] pts, vec[] nAra, vec[] tAra, boolean derived) {
//        for(int i = 0; i < pts.length; ++i){
//            pts[i].drawNorm((derived ? 0 : 1), nAra[i], tAra[i]);
//        }
//    }
    public void drawCntlPts(cntlPt[] pts, vec[] nAra, vec[] tAra, boolean derived) {
    	pa.pushStyle();
    	int clrInt = 0;
        for(int i = 0; i < pts.length; ++i){
        	clrInt = (int)(i/(1.0f * pts.length) * 255.0f);
            pa.fill(clrInt,0,(255 - clrInt),255);  
            pa.stroke(clrInt,0,(255 - clrInt),255); 
            pts[i].drawRad(nAra[i], tAra[i]);
        }
        pa.popStyle();
    }
        



}//myNormOffset


