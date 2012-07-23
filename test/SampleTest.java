import com.affymetrix.genometryImpl.AnnotatedSeqGroup;
import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.symloader.BED;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import findjunction.FindJunction;
import java.io.File;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 *
 * @author hiralv
 */
public class SampleTest {
    
    String input,output;
    
    @Before
    public void setUp() throws URISyntaxException, Exception{
        input = "data/small_hits.bam";
        output ="data/small_hits.bed";
        FindJunction fJ = new FindJunction();
        fJ.init(input, output, 5, false, null, true);
    }
    
    @Test
    public void sample() throws URISyntaxException, Exception{
        File outputFile = new File(output);
        outputFile = outputFile.getAbsoluteFile();
        output = outputFile.getAbsolutePath();
        URI bedUri = new URI("file:"+output);
        Map<String, SeqSymmetry> map = new HashMap<String, SeqSymmetry>();
        AnnotatedSeqGroup group = new AnnotatedSeqGroup("small_hits");
        BED bedFile = new BED(bedUri,"small_hits",group);
        List<BioSeq> chromosomes = bedFile.getChromosomeList();
        for(BioSeq chr : chromosomes){
            List<SeqSymmetry> symList = bedFile.getChromosome(chr); 
            for(SeqSymmetry sym : symList){
                SeqSpan span = sym.getSpan(chr);
                String key = ""+span.getMin()+""+span.getMax();
                Assert.assertEquals(false, map.containsKey(key));
                map.put(key, sym);
            }
            map.clear();
            symList.clear();
        }
    }
    
    @After
    public void end(){
       File file = new File(output);
       file.delete();
    }
}
