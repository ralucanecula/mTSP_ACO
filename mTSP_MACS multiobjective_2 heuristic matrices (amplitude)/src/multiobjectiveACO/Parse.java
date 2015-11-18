package multiobjectiveACO;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class Parse {

    static class OptComparator implements Comparator<Option> {

		Map<String, Integer> opt = new HashMap<String, Integer>();
	
		public OptComparator() {
		    int i = 0;
		    
		    opt.put("u", i++);
		    opt.put("z", i++);
		}
	
		@Override
		public int compare(Option o1, Option o2) {
		    if (o1.getValue() == null || o2.getValue() == null) return 0;
		    else
		    	return (opt.get(o1.getOpt()) - opt.get(o2.getOpt()));
		}
    }

    static void parse_commandline(String args[]) {
		if (args.length == 0) {
		    System.err.println("No options are specified.");
		    System.err.println("Try `--help' for more information.");
		    System.exit(1);
		}
	
		Options options = new Options();
		options.addOption("u", "as", false, "apply basic Ant System");
		options.addOption("z", "acs", false, "apply ant colony colony system");
	
		CommandLine cmd = null;
		CommandLineParser parser = new BasicParser();
		try {
		    cmd = parser.parse(options, args);
		} catch (ParseException e) {
		    System.err.println("Error: " + e.getMessage());
		    System.exit(1);
		}
	
		// Choice of ONE algorithm
		int algorithmCount = 0;
		if (cmd.hasOption("u")) {
		    algorithmCount++;
		}
		if (cmd.hasOption("z")) {
		    algorithmCount++;
		}
		if (algorithmCount > 1) {
		    System.err.println("Error: More than one ACO algorithm enabled in the command line.");
		    System.exit(1);
		} else if (algorithmCount == 1) {
		    Ants.as_flag = false;
		    Ants.acs_flag = false;
		}
	
		if (cmd.hasOption("u")) {
		    Ants.as_flag = true;
		    InOut.set_default_as_parameters();
		    System.out.println("\nRun basic Ant System");
		}
		if (cmd.hasOption("z")) {
		    Ants.acs_flag = true;
		    InOut.set_default_acs_parameters();
		    System.out.println("\nRun Ant Colony System");
		}

    }
}
