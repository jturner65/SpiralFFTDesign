package SpiralFFTDesPKG;

import SpiralFFTDesPKG.FFTSpiralDsn.*;
//use to denote control points for drawn strokes
public class cntlPt extends pt {
		public static FFTSpiralDsn pa;
		public int ID;
		public static int IDincr = 0;
		public static final float maxR = 75, 
				minR = 1,
				baseRad = 20;			//default radius for control points
		public float r, w;				//weight is calculated based on the distance to neighboring cntl pts when cntl pts are drawn
		public static int[][] clr = new int[][]{{0,140,240,255}, {0,0,0,255}};
		
		public cntlPt(FFTSpiralDsn _pa, pt _p, float _r, float _w){ _pa.super(_p.x,_p.y); pa=_pa; ID=IDincr++;r=_r; w=_w; }
		public cntlPt(FFTSpiralDsn _pa, pt _p, float _w){this(_pa, _p, baseRad, _w);}
		public cntlPt(FFTSpiralDsn _pa, pt _p){this(_pa, _p, baseRad, baseRad);}
		public cntlPt(FFTSpiralDsn _pa){this(_pa, _pa.P(0,0),1);}
		public cntlPt(cntlPt _p){this(cntlPt.pa, cntlPt.pa.P(_p),_p.w); r = _p.r; w = _p.w;ID = _p.ID;}		
		public static cntlPt L(cntlPt A, float s, cntlPt B){	return new cntlPt(cntlPt.pa, cntlPt.pa.L((pt)A, s, (pt)B), capInterpR(A.r, s, B.r), (1-s)*A.w + (s)*B.w);}//(1-s)*A.r + (s)*B.r,
		public static cntlPt P(cntlPt A, cntlPt B){	float s = .5f;return L(A, s, B);}
		public pt set(pt P){super.setTo(P); return (pt)this;}
		private static float capInterpR(float a, float s, float b){ float res = (1-s)*a + (s)*b; res = (res < minR ? minR : res > maxR ? maxR : res); return res;}
		public void drawMe(int cIdx){	pa.fill(clr[cIdx][0],clr[cIdx][1],clr[cIdx][2],clr[cIdx][3]);  pa.stroke(clr[cIdx][0],clr[cIdx][1],clr[cIdx][2],clr[cIdx][3]);		super.show(2);}		
		public void drawRad(int cIdx,vec I, vec J){
            pa.fill(clr[cIdx][0],clr[cIdx][1],clr[cIdx][2],clr[cIdx][3]);  
            pa.stroke(clr[cIdx][0],clr[cIdx][1],clr[cIdx][2],clr[cIdx][3]); 
            pa.circle(this, r, I,J,20);
	    }
		public void drawRad(vec I, vec J){
            pa.circle(this, r, I,J,20);
	    }
		public void drawNorm(int cIdx,vec I, vec J) {
		    pt p1 = pa.P(this).add(pa.V(I).mul(r));
            pt p2 = pa.P(this).add(pa.V(I).mul(-r));
            pa.stroke(0,0,0,255);
            pa.edge(p1, p2); 
		}
		public void calcRadFromWeight(float lenRatio, boolean inv){r = Math.min(maxR, Math.max(minR, baseRad * (inv ?  (lenRatio/w) : (pa.wScale*w/(lenRatio*lenRatio)))));  }
		public void modRad(float modAmt){float oldR = r; r += modAmt; r = (r < minR ? minR : r > maxR ? maxR : r); w *= oldR/r; }
		public String toString(){String res = "Cntl Pt ID:"+ID+" p:"+super.toString()+" r:"+r+" w:"+w;return res;}
}//class cntlPoint

