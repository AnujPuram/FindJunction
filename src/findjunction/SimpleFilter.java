/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction;

import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.filter.SymmetryFilterI;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;

/**
 *
 * @author auser
 */
public class SimpleFilter implements SymmetryFilterI{

    @Override
    public String getName() {
        return "simple";
    }

    @Override
    public boolean setParam(Object o) {
        return true;
    }

    @Override
    public Object getParam() {
        return null;
    }

    @Override
    public boolean filterSymmetry(BioSeq bioseq, SeqSymmetry ss) {
        SeqSpan span = ss.getSpan(bioseq);
        if(Math.abs(span.getEnd()-span.getStart()) >= 100)
            return true;
        return false;
    }
    
}
