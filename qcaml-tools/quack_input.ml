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
  (*
  let open Printf in
  *)
  let g_basis = Qcaml.Gaussian.Basis.general_basis basis in
  let dict = 
    Array.to_list nuclei
    |> List.mapi (fun i (e, _) ->
          (i+1), List.assoc e g_basis
        ) 
  in
  List.iter (fun (i,b) ->
     Format.printf "%d  %d\n" i (Array.length b);
     Array.iter (fun x ->
        Format.printf "%a" Qcaml.Gaussian.General_basis.pp_gcs x) b
     ) dict
(*
  List.iteri (fun atom_number (_element, _basis) ->
    printf "%3d %3d\n" (atom_number+1) 0;
    ) basis
    *)
    




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
    Qcaml.Particles.Nuclei.of_xyz_file nuclei_file
  in

  let basis = 
    Qcaml.Gaussian.Basis.of_nuclei_and_basis_filename ~nuclei basis_file
  in

  (* Print basis *)
  print_basis nuclei basis;
  ()

let () =
  let usage_msg = "Available options:" in
  Arg.parse speclist (fun _ -> ()) usage_msg;
  run ()


