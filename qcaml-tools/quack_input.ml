let quack_dir = 
  try Sys.getenv "QUACK_DIR" with
  Not_found -> "."

let quack_basis_filename = quack_dir ^ "/input/basis"
let quack_molecule_filename = quack_dir ^ "/input/molecule"

module Command_line = Qcaml.Common.Command_line
module Util = Qcaml.Common.Util

let () =
  let open Command_line in
  begin
    set_header_doc (Sys.argv.(0) ^ " - QuAcK command");
    set_description_doc "Prepares the input data for QuAcK.
Writes $QUACK_DIR/input/basis and $QUACK_DIR/input/molecule.
If $QUACK_DIR is not set, $QUACK_DIR is replaces by the current
directory.";
    set_specs
      [ { short='b' ; long="basis" ; opt=Mandatory;
          arg=With_arg "<string>";
          doc="Name of the file containing the basis set"; } ;

        { short='x' ; long="xyz" ; opt=Mandatory;
          arg=With_arg "<string>";
          doc="Name of the file containing the nuclear coordinates in xyz format"; } ;

        { short='m' ; long="multiplicity" ; opt=Optional;
          arg=With_arg "<int>";
          doc="Spin multiplicity (2S+1). Default is singlet"; } ;

        { short='c' ; long="charge" ; opt=Optional;
          arg=With_arg "<int>";
          doc="Total charge of the molecule. Specify negative charges with 'm' instead of the minus sign, for example m1 instead of -1. Default is 0"; } ;

        { short='f' ; long="frozen-core" ; opt=Optional;
          arg=Without_arg ;
          doc="Freeze core MOs. Default is false"; } ;

        { short='r' ; long="rydberg" ; opt=Optional;
          arg=With_arg "<int>" ;
          doc="Number of Rydberg electrons. Default is 0"; } ;
      ]
  end;

  (* Handle options *)
  let basis_file  = Util.of_some @@ Command_line.get "basis" in
  let nuclei_file = Util.of_some @@ Command_line.get "xyz" in
  let frozen_core = Command_line.get_bool "frozen-core" in

  let charge =
    match Command_line.get "charge" with
    | Some x -> ( if x.[0] = 'm' then
                    ~- (int_of_string (String.sub x 1 (String.length x - 1)))
                  else
                    int_of_string x )
    | None   -> 0
  in

  let multiplicity =
    match Command_line.get "multiplicity" with
    | Some x -> int_of_string x
    | None -> 1
  in

  let rydberg =
    match Command_line.get "rydberg" with
    | Some x -> int_of_string x
    | None -> 0
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
  in
  print_basis nuclei basis;


  let print_molecule nuclei electrons =
    let oc = open_out quack_molecule_filename in
    let nat  = Array.length nuclei in
    let nela = Qcaml.Particles.Electrons.n_alfa electrons in
    let nelb = Qcaml.Particles.Electrons.n_beta electrons in
    let ncore = 
      if frozen_core then
        Qcaml.Particles.Nuclei.small_core nuclei 
      else 0
    in
    let nryd = rydberg in
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
  in    
  print_molecule nuclei electrons;

  ()



