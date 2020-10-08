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

let run () =
  let basis_file =
    match !basis_file with
    | None -> raise (Invalid_argument "Basis set file should be specified with -b")
    | Some x -> x
  and nuclei_file =
    match !nuclei_file with
    | None -> raise (Invalid_argument "Coordinate file should be specified with -x")
    | Some x -> x
  in

  let nuclei = 
    Qcaml.Nuclei.of_xyz_file nuclei_file
  in

  let basis = 
    QCaml.Gaussian_basis.of_nuclei_and_basis_filename  nuclei  basis_file
    |> QCaml.Gaussian_basis.general_basis
  in

  (* Print basis *)
  Format.printf "%a" QCaml.Gaussian_basis.pp basis;
  ()
(*
  List.map (fun (element, shell) ->
    Simulation.of_filenames ?range_separation ?charge ?multiplicity
      ~nuclei:nuclei_file basis_file
  in

  print_endline @@ Nuclei.to_string @@ Simulation.nuclei s;
  print_endline "Nuclear repulsion : ";
  print_float @@ Simulation.nuclear_repulsion s; print_newline ();
  print_endline @@ Basis.to_string  @@ Simulation.basis s;

  let ao_basis = Simulation.ao_basis s in
  let overlap  = AOBasis.overlap  ao_basis in
  let eN_ints  = AOBasis.eN_ints  ao_basis in
  let kin_ints = AOBasis.kin_ints ao_basis in
  let ee_ints  = AOBasis.ee_ints  ao_basis in
  Overlap.to_file ~filename:("Ov.dat") overlap;
  NucInt.to_file ~filename:("Nuc.dat") eN_ints;
  KinInt.to_file ~filename:("Kin.dat") kin_ints;
  ERI.to_file    ~filename:("ERI.dat") ee_ints;
  match range_separation with
  | Some _mu ->
      ERI_lr.to_file ~filename:("ERI_lr.dat") (AOBasis.ee_lr_ints ao_basis)
  | None -> ()
*)

let () =
  let usage_msg = "Available options:" in
  Arg.parse speclist (fun _ -> ()) usage_msg;
  run ()


