/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction.filters;

import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.filter.SymmetryFilterI;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import findjunction.FindJunction;

/**
 *
 * @author auser
 */
public class DuplicateSymFilter implements SymmetryFilterI{

    private static int start;
    @Override
    public String getName() {
        return null;
    }

    @Override
    public boolean setParam(Object o) {
        start = (Integer)o;
        return true;
    }

    @Override
    public Object getParam() {
        return start;
    }

    @Override
    public boolean filterSymmetry(BioSeq bioseq, SeqSymmetry ss) {
        if(ss.getSpan(bioseq).getMin() >= (Integer)getParam())
            return true;
        else
            return false;
    }
    
}
