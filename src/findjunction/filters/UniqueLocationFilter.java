/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction.filters;

import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.filter.SymmetryFilterI;
import com.affymetrix.genometryImpl.symmetry.BAMSym;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import com.affymetrix.genometryImpl.symmetry.SymWithProps;

/**
 *
 * @author auser
 */
public class UniqueLocationFilter implements SymmetryFilterI{

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
    public boolean filterSymmetry(BioSeq bioseq, SeqSymmetry ss) {
        int currentNH = (Integer)(((BAMSym)ss).getProperty("NH"));
        if(currentNH != 1)
            return false;
        return true;
    }
    
}
