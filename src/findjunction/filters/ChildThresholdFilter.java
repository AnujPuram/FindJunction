/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction.filters;

import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.filter.SymmetryFilterI;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import findjunction.FindJunction;

/**
 *
 * @author Anuj
 */
public class ChildThresholdFilter implements SymmetryFilterI{

    int threshold;
    
    @Override
    public String getName() {
        return null;
    }

    @Override
    public boolean setParam(Object o) {
        threshold = (Integer)o;
        return true;
    }

    @Override
    public Object getParam() {
        return threshold;
    }

    @Override
    public boolean filterSymmetry(BioSeq bioseq, SeqSymmetry ss) {
        SeqSpan span = ss.getSpan(bioseq);
        if((span.getMax() - span.getMin()) < threshold)
            return false;
        else
            return true;
    }
}
