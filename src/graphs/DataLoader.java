package graphs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class DataLoader {
	private static String data_file = "./data/enron-mysqldump.sql";
	public static final int ADVOGATO 		  = 0;
	public static final int ENRON_MULTI_EDGE  = 1;
	public static final int ENRON_SINGLE_EDGE = 2;
	public static final int CONGRESS_TWITTER  = 3;
	
	private static HashMap<String, Integer> data_set_ids_by_name = new HashMap<String, Integer>(10);
	static{
		data_set_ids_by_name.put("ADVOGATO", ADVOGATO);
		data_set_ids_by_name.put("ADVOGATO".toLowerCase(), ADVOGATO);
		data_set_ids_by_name.put("ENRON", ENRON_SINGLE_EDGE);
		data_set_ids_by_name.put("ENRON".toLowerCase(), ENRON_SINGLE_EDGE);
		data_set_ids_by_name.put("CONGRESS", CONGRESS_TWITTER);
		data_set_ids_by_name.put("CONGRESS".toLowerCase(), CONGRESS_TWITTER);
	}
	
	public static Graph get_graph(String name){
		return get_graph(data_set_ids_by_name.get(name));
	}
	
	public static Graph get_graph(Integer id){
		if(id==ADVOGATO) {
			return get_advogato_graph();
		}else if(id==ENRON_MULTI_EDGE) {
			return get_enron_graph();
		}else if(id==ENRON_SINGLE_EDGE) {
			Graph g = DataLoader.get_enron_single_edge_graph();
			return g;
		}else if(id==CONGRESS_TWITTER) {
			return get_congress_twitter_graph();
		}else{
			System.err.println("get_graph("+id+") Unkown id");
			return null;
		}
	}
	
	public static String get_graph_name(int graph_id) {
		if(graph_id==CONGRESS_TWITTER) {
			return "CONGRESS_TWITTER";
		}else if(graph_id==ADVOGATO){
			return "ADVOGATO";
		}else if(graph_id==ENRON_SINGLE_EDGE){
			return "ENRON";
		}else{
			return "Unknown graph id="+graph_id;
		}
	}
	
	public static String get_synthetic_path(String graph_name) {
		return "./data/synthetic/"+graph_name;
	}
	
	public static ArrayList<MaterializedGraphs> get_sanitized_graphs(Integer id){
		String graph_name = get_graph_name(id);
		String path = get_synthetic_path(graph_name);
		Set<String> all_files = null;
		try {
			all_files = listFilesUsingFilesList(path);
			ArrayList<MaterializedGraphs> mg = new ArrayList<MaterializedGraphs>(all_files.size());
			for(String s : all_files) {
				mg.add(new MaterializedGraphs(s, graph_name));
			}
			Collections.sort(mg);
			return mg;
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	public static Set<String> listFilesUsingFilesList(String dir) throws IOException {
	    try (Stream<Path> stream = Files.list(Paths.get(dir))) {
	        return stream
	          .filter(file -> !Files.isDirectory(file))
	          .map(Path::getFileName)
	          .map(Path::toString)
	          .collect(Collectors.toSet());
	    }
	}
	
	public static String get_graph_name(Integer id){
		if(id==ADVOGATO) {
			return "ADVOGATO";
		}else if(id==ENRON_MULTI_EDGE) {
			return "ENRON_MULTI_EDGE".toLowerCase();
		}else if(id==ENRON_SINGLE_EDGE) {
			return "ENRON".toLowerCase();
		}else if(id==CONGRESS_TWITTER) {
			return "CONGRESS_TWITTER".toLowerCase();
		}else{
			System.err.println("get_graph("+id+") Unkown id");
			return null;
		}
	}
	
	static Graph get_enron_graph(){
		final String is_message = "INSERT INTO message VALUES (";
		final String is_recipient = "INSERT INTO recipientinfo VALUES (";
		final int message_mid 		= 0;
		final int message_sender 	= 1;
		final int recipient_mid		= 1;
		final int recipient_email 	= 3;
		
		/**
		 * [mid,email]
		 */
		ArrayList<String[]> from = new ArrayList<String[]>(100000);
		HashMap<String, ArrayList<String>> to = new HashMap<String, ArrayList<String>>(100000);
		
		File f = new File(data_file);
		if(f.exists()) {
			try(BufferedReader br = new BufferedReader(new FileReader(f))){
				String line;
				int i = 0;
				while((line = br.readLine()) != null) {
					if(line.startsWith(is_message)) {
						String[] tokens = line.replace(is_message, "").split(",");
						String[] temp = {tokens[message_mid],tokens[message_sender]};
						from.add(temp);
						//System.out.println(tokens[0]+" "+tokens[1]);
					}else if(line.startsWith(is_recipient)){
						String[] tokens = line.replace(is_recipient, "").split(",");
						ArrayList<String> all_recipients = to.get(tokens[recipient_mid]);
						if(all_recipients==null) {
							all_recipients = new ArrayList<String>();
							to.put(tokens[recipient_mid], all_recipients);
						}
						all_recipients.add(tokens[recipient_email]);
						//System.out.println(tokens[1]+" "+tokens[3]);
					}
					
					if(i++%100000==0)	
						System.out.println(i+" "+line);
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			// get all nodes
			HashSet<String> all_email_adresses = new HashSet<String>();
			for(String[] message : from) {
				all_email_adresses.add(message[1]);
			}
			for( Entry<String, ArrayList<String>> entry : to.entrySet()) {
				for(String email : entry.getValue()) {
					all_email_adresses.add(email);	
				}
			}
			System.out.println("Got |G| = "+all_email_adresses.size());
			String[] nodes = new String[all_email_adresses.size()];
			int counter = 0;
			for(String s : all_email_adresses) {
				nodes[counter++] = s;
			}
			Arrays.sort(nodes);
			HashMap<String, Integer> node_to_id = new HashMap<String, Integer>(nodes.length);
			for(int id=0;id<nodes.length;id++) {
				node_to_id.put(nodes[id], id);
			}
			
			Graph g = new Graph(all_email_adresses.size(), "enron");
			for(String[] message : from) {
				String mid = message[0];
				final int start_vertex = node_to_id.get(message[1]);
				ArrayList<String> recipients = to.get(mid);
				if(recipients!=null) {
					for(String node : recipients) {
						int target_vertex = node_to_id.get(node);
						//System.out.println(start_vertex+" "+target_vertex);
						g.add_edge(start_vertex, target_vertex);
					}
				}else{
					//System.err.println("recipients==null "+Arrays.toString(message));
				}
			}
			System.out.println(g);
			return g;
		}else{
			System.err.println(f+" does not exist");
			return null;
		}
	}
	
	/**
	 * 
	 * asserts that employees come first and messages second.
	 * @return
	 */
	@Deprecated
	static Graph get_enron_graph_employee_unique_edges(){
		final String is_employee = "INSERT INTO employeelist VALUES (";
		final String is_message = "INSERT INTO message VALUES (";
		final String is_recipient = "INSERT INTO recipientinfo VALUES (";
		final int message_mid 		= 0;
		final int message_sender 	= 1;
		final int recipient_mid		= 1;
		final int recipient_email 	= 3;
		final int employee_email_id = 3;
		
		/**
		 * [mid,email]
		 */
		ArrayList<String[]> from = new ArrayList<String[]>(100000);
		HashMap<String, ArrayList<String>> to = new HashMap<String, ArrayList<String>>(100000);
		
		File f = new File(data_file);
		if(f.exists()) {
			try(BufferedReader br = new BufferedReader(new FileReader(f))){
				String line;
				int i = 0;
				/**
				 * Only of employees
				 */
				HashSet<String> email_ids = new HashSet<String>();
				HashSet<String> mids = new HashSet<String>();
				
				while((line = br.readLine()) != null) {
					if(line.startsWith(is_message)) {
						String[] tokens = line.replace(is_message, "").split(",");
						if(email_ids.contains(tokens[message_sender])) {//only messages from enron employees
							String[] temp = {tokens[message_mid],tokens[message_sender]};
							from.add(temp);	
							mids.add(tokens[message_mid]);
						}
						//System.out.println(tokens[0]+" "+tokens[1]);
					}else if(line.startsWith(is_recipient)){
						String[] tokens = line.replace(is_recipient, "").split(",");
						if(mids.contains(tokens[recipient_mid])) {//Mail from an enron employee
							if(email_ids.contains(tokens[recipient_mid])) {//Mail to an enron employee
								ArrayList<String> all_recipients = to.get(tokens[recipient_mid]);
								if(all_recipients==null) {
									all_recipients = new ArrayList<String>();
									to.put(tokens[recipient_mid], all_recipients);
								}
								all_recipients.add(tokens[recipient_email]);
								//System.out.println(tokens[1]+" "+tokens[3]);
							}
						}
					}else if(line.startsWith(is_employee)){
						String[] tokens = line.replace(is_employee, "").split(",");
						String email_id = tokens[employee_email_id].replace(");", "");
						//System.out.println(email_id);
						email_ids.add(email_id);
					}
					
					if(i++%100000==0)	
						System.out.println(i+" "+line);
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			// get all nodes
			HashSet<String> all_email_adresses = new HashSet<String>();
			for(String[] message : from) {
				all_email_adresses.add(message[message_sender]);
			}
			for( Entry<String, ArrayList<String>> entry : to.entrySet()) {
				for(String email : entry.getValue()) {
					all_email_adresses.add(email);	
				}
			}
			System.out.println("Got |G| = "+all_email_adresses.size());
			String[] nodes = new String[all_email_adresses.size()];
			int counter = 0;
			for(String s : all_email_adresses) {
				nodes[counter++] = s;
			}
			Arrays.sort(nodes);
			HashMap<String, Integer> node_to_id = new HashMap<String, Integer>(nodes.length);
			for(int id=0;id<nodes.length;id++) {
				node_to_id.put(nodes[id], id);
			}
			
			Graph g = new Graph(all_email_adresses.size(),"enron unique e");
			for(String[] message : from) {
				String mid = message[0];
				final int start_vertex = node_to_id.get(message[1]);
				ArrayList<String> recipients = to.get(mid);
				if(recipients!=null) {
					for(String node : recipients) {
						int target_vertex = node_to_id.get(node);
						//System.out.println(start_vertex+" "+target_vertex);
						g.add_edge(start_vertex, target_vertex);
					}
				}else{
					//System.err.println("recipients==null "+Arrays.toString(message));
				}
			}
			System.out.println(g);
			return g;
		}else{
			System.err.println(f+" does not exist");
			return null;
		}
	}
	
	private static String advogato_edges = "./data/advogato_edges.csv";
	private static String advogato_nodes = "./data/advogato_nodes.csv";
	
	static Graph get_advogato_graph() {
		final int source_offset = 0;
		final int target_offset = 1;
		final int num_node = 6541;//According to advogato_nodes [0,6540]
		Graph g = new Graph(num_node, "advogato");
		
		File f = new File(advogato_edges);
		if(f.exists()) {
			try(BufferedReader br = new BufferedReader(new FileReader(f))){
				String line;
				br.readLine();//skip meta info
				
				while((line = br.readLine()) != null) {
					String tokens[] = line.split(",");
					int source = Integer.parseInt(tokens[source_offset]);
					int target = Integer.parseInt(tokens[target_offset]);
					if(source<num_node && source>=0) {
						if(target<num_node && target>=0) {
							g.add_edge(source, target);
						}else{
							System.err.println("target<num_node && target>=0");
						}
					}else{
						System.err.println("source<num_node && source>=0");
					}
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			System.out.println(g);
			return g;
		}else{
			System.err.println(f+" does not exist");
			return null;
		}
	}
	
	static Graph get_graph_from_edge_list(String file, int num_node) {
		final int source_offset = 0;
		final int target_offset = 1;
		Graph g = new Graph(num_node, file);
		
		File f = new File(file);
		if(f.exists()) {
			try(BufferedReader br = new BufferedReader(new FileReader(f))){
				String line;
				
				while((line = br.readLine()) != null) {
					try {
						String tokens[] = line.split(" ");
						int source = Integer.parseInt(tokens[source_offset]);
						int target = Integer.parseInt(tokens[target_offset]);
						if(source<num_node && source>=0) {
							if(target<num_node && target>=0) {
								g.add_edge(source, target);
							}else{
								System.err.println("get_graph_from_edge_list() target<num_node && target>=0");
							}
						}else{
							System.err.println("get_graph_from_edge_list() source<num_node && source>=0");
						}
					} catch (Exception e) {
						// TODO: handle exception
					}
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			System.out.println(g);
			return g;
		}else{
			System.err.println(f+" does not exist");
			return null;
		}
	}
	private static String congress_twitter = "./data/congress.edgelist";
	static Graph get_congress_twitter_graph() {
		final int num_node = 475;//According to https://snap.stanford.edu/data/congress-twitter.html
		return get_graph_from_edge_list(congress_twitter, num_node);
	}
	private static String ennon_single_edge = "./data/dp_enron_dedup.edgelist";
	static Graph get_enron_single_edge_graph() {
		final int num_node = 75557;
		return get_graph_from_edge_list(ennon_single_edge, num_node);
	}
	
	public static Graph load(MaterializedGraphs mg) {
		int num_node = -1;
		if(mg.graph_name.equals(get_graph_name(0))) {
			num_node = 6541;//Advogato
		}else if(mg.graph_name.equals(get_graph_name(2).toLowerCase())) {
			num_node = 75557;//Enron
		}else if(mg.graph_name.equals(get_graph_name(3))) {
			num_node = 475;//congress_twitter According to https://snap.stanford.edu/data/congress-twitter.html
		}else if(mg.graph_name.equals(get_graph_name(3).toLowerCase())) {
			num_node = 475;//congress_twitter According to https://snap.stanford.edu/data/congress-twitter.html
		}else{
			
		}
		String path = get_synthetic_path(mg.graph_name)+"/"+mg.path;
		return get_graph_from_edge_list(path, num_node);
	}

	public static ArrayList<Graph> get_sanitized_graphs(ArrayList<MaterializedGraphs> mg_s) {
		ArrayList<Graph> ret = new ArrayList<Graph>(mg_s.size());
		for(MaterializedGraphs mg : mg_s) {
			ret.add(load(mg));
		}
		return ret;
	}
}
