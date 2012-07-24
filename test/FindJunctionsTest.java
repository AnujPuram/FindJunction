import com.affymetrix.genometryImpl.AnnotatedSeqGroup;
import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.symloader.BED;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import findjunction.FindJunction;
import java.io.File;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
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
public class FindJunctionsTest {
    
    String input;
    List<String> output = new ArrayList<String>();
    
    @Before
    public void setUp() throws URISyntaxException, Exception{
        File inputDirectory = new File("data/input");
        inputDirectory = inputDirectory.getAbsoluteFile();
        FindJunction fJ = new FindJunction();
        if(inputDirectory.isDirectory()){
            File files[] = inputDirectory.listFiles();
            for(int i=0;i<files.length;i++){
                String name = files[i].getName();
                String ext = name.substring(name.length()-3, name.length());
                if(ext.equals("bam")){
                    String out = name.substring(0, name.length()-3)+"bed";
                    fJ.init("data/input/"+name, "data/output/"+out, 5, false, null, true);
                    output.add(out);
                }
            }
        }
        
    }
    
    @Test
    public void sample() throws URISyntaxException, Exception{
        File outputDirectory = new File("data/output");
        outputDirectory = outputDirectory.getAbsoluteFile();
        if(outputDirectory.isDirectory()){
            File outputFiles[] = outputDirectory.listFiles();
            for(File file : outputFiles){
                String name = file.getName();
                String ext = name.substring(name.length()-3, name.length());
                if(!ext.equals("bed")){
                    continue;
                }
                file = file.getAbsoluteFile();
                String path = file.getAbsolutePath();
                URI bedUri = new URI("file:"+path);
                Map<String, SeqSymmetry> map = new HashMap<String, SeqSymmetry>();
                AnnotatedSeqGroup group = new AnnotatedSeqGroup(name.substring(0 , name.length()-3));
                BED bedFile = new BED(bedUri,name.substring(0 , name.length()-3),group);
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
                System.out.println(file.getName()+" success....");
            }
        }
    }
    
    @After
    public void end(){
       File outputDirectory = new File("data/output");
       outputDirectory = outputDirectory.getAbsoluteFile();
       if(outputDirectory.isDirectory()){
           for(File file : outputDirectory.listFiles()){
               if(output.contains(file.getName())){
                   file.delete();
               }
           }
       }
       output.clear();
    }
}