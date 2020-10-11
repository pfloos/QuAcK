let out_file         : string option ref = ref None
let basis_file       : string option ref = ref None
let nuclei_file      : string option ref = ref None
let charge           : int    option ref = ref None
let multiplicity     : int    option ref = ref None
let range_separation : float  option ref = ref None


let speclist = [
  ( "-b" , Arg.String (fun x -> basis_file := Some x),
    "File containing the atomic basis set") ;
  ( "-x" , Arg.String (fun x -> nuclei_file := Some x),
    "File containing the nuclear coordinates") ;
  ( "-u" , Arg.Float  (fun x -> range_separation := Some x),
    "Value of mu, the range separation factor") ;
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
  and range_separation = !range_separation
  in

  let nuclei =
    Qcaml.Particles.Nuclei.of_xyz_file nuclei_file
  in

  let operators = 
    match range_separation with
    | None -> []
    | Some mu -> [ Qcaml.Operators.Operator.of_range_separation mu ]
  in

  let ao_basis =
    Qcaml.Ao.Basis.of_nuclei_and_basis_filename ~kind:`Gaussian
      ~operators ~cartesian:true ~nuclei basis_file
  in

  let overlap  = Qcaml.Ao.Basis.overlap  ao_basis in
  let eN_ints  = Qcaml.Ao.Basis.eN_ints  ao_basis in
  let kin_ints = Qcaml.Ao.Basis.kin_ints ao_basis in
  let ee_ints  = Qcaml.Ao.Basis.ee_ints  ao_basis in
  Qcaml.Gaussian_integrals.Overlap.to_file ~filename:("Ov.dat") overlap;
  Qcaml.Gaussian_integrals.Electron_nucleus.to_file ~filename:("Nuc.dat") eN_ints;
  Qcaml.Gaussian_integrals.Kinetic.to_file ~filename:("Kin.dat") kin_ints;
  Qcaml.Gaussian_integrals.Eri.to_file    ~filename:("ERI.dat") ee_ints;
  match range_separation with
  | Some _mu ->
      Qcaml.Gaussian_integrals.Eri_long_range.to_file ~filename:("ERI_lr.dat") (Qcaml.Ao.Basis.ee_lr_ints ao_basis)
  | None -> ()


let () =
  let usage_msg = "Available options:" in
  Arg.parse speclist (fun _ -> ()) usage_msg;
  run ()

