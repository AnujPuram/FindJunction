/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction.filters;

import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.filter.SymmetryFilterI;
import com.affymetrix.genometryImpl.symmetry.BAMSym;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import findjunction.FindJunction;
import net.sf.samtools.Cigar;

/**
 *
 * @author Anuj
 */
public class ThresholdFilter implements SymmetryFilterI{

    public FindJunction findJunction = new FindJunction();
    
    @Override
    public String getName() {
        return null;
    }

    @Override
    public boolean setParam(Object o) {
        return false;
    }

    @Override
    public Object getParam() {
        return null;
    }

    @Override
    //This filter is going to discard the whole sym even if one of its children has less length than the threshold
    public boolean filterSymmetry(BioSeq bioseq, SeqSymmetry ss){
        int threshold = findJunction.getThreshold(); 
        int childCount = ss.getChildCount();
        ChildThresholdFilter filter = new ChildThresholdFilter();
        SeqSymmetry child;
        for(int i=0;i<childCount;i++){
            child = ss.getChild(i);
            if(!filter.filterSymmetry(bioseq, child))
                return false;
        }
        return true;
    }   
}