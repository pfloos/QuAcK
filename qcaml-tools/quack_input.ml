let quack_dir = 
  try Sys.getenv "QUACK_DIR" with
  Not_found -> "."

let quack_basis_filename = quack_dir ^ "/input/basis"
let quack_molecule_filename = quack_dir ^ "/input/molecule"


let basis_file       : string option ref = ref None
let nuclei_file      : string option ref = ref None
let charge           : int    option ref = ref None
let multiplicity     : int    option ref = ref None


let speclist = [
  ( "-b" , Arg.String (fun x -> basis_file := Some x),
    "File containing the atomic basis set") ;
  ( "-c" , Arg.Int    (fun x -> charge := Some x),
    "Total charge of the system") ;
  ( "-m" , Arg.Int    (fun x -> multiplicity := Some x),
    "Multiplicity of the system") ;
  ( "-x" , Arg.String (fun x -> nuclei_file := Some x),
    "File containing the nuclear coordinates") ;
]


let print_basis nuclei basis =
  let oc = open_out quack_basis_filename in
  let ocf = Format.formatter_of_out_channel oc in
  let g_basis = Qcaml.Gaussian.Basis.general_basis basis in
  let dict = 
    Array.to_list nuclei
    |> List.mapi (fun i (e, _) ->
          (i+1), List.assoc e g_basis
        ) 
  in
  List.iter (fun (i,b) ->
     Format.fprintf ocf "%d  %d\n" i (Array.length b);
     Array.iter (fun x ->
        Format.fprintf ocf "%a" Qcaml.Gaussian.General_basis.pp_gcs x) b
     ) dict;
  close_out oc
    


let print_molecule nuclei electrons =
  let oc = open_out quack_molecule_filename in
  let nat  = Array.length nuclei in
  let nela = Qcaml.Particles.Electrons.n_alfa electrons in
  let nelb = Qcaml.Particles.Electrons.n_beta electrons in
  let ncore = Qcaml.Particles.Nuclei.small_core nuclei in
  let nryd  = 0 in
  Printf.fprintf oc "# nAt nEla nElb nCore nRyd\n";
  Printf.fprintf oc " %4d %4d %4d %5d %4d\n" nat nela nelb ncore nryd;
  Printf.fprintf oc "# Znuc   x            y           z\n";
  let open Qcaml.Common.Coordinate in
  Array.iter (fun (element, coord) ->
    Printf.fprintf oc "%3s    %16.10f     %16.10f     %16.10f\n"
        (Qcaml.Particles.Element.to_string element)
        coord.x coord.y coord.z
      ) nuclei;
  close_out oc
    



let run () =
  let basis_file =
    match !basis_file with
    | None -> raise (Invalid_argument "Basis set file should be specified with -b")
    | Some x -> x
  and nuclei_file =
    match !nuclei_file with
    | None -> raise (Invalid_argument "Coordinate file should be specified with -x")
    | Some x -> x
  and charge = 
    match !charge with
    | None -> 0
    | Some c -> c
  and multiplicity = 
    match !multiplicity with
    | None -> 1
    | Some m -> m
  in

  let nuclei = 
    Qcaml.Particles.Nuclei.of_xyz_file nuclei_file
  in

  let electrons =
    Qcaml.Particles.Electrons.of_atoms ~multiplicity ~charge nuclei
  in

  let basis = 
    Qcaml.Gaussian.Basis.of_nuclei_and_basis_filename ~nuclei basis_file
  in

  (* Print basis *)
  print_molecule nuclei electrons;
  print_basis nuclei basis;
  ()

let () =
  let usage_msg = "Available options:" in
  Arg.parse speclist (fun _ -> ()) usage_msg;
  run ()


