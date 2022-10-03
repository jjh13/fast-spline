#!/usr/bin/env python

import sys
import argparse
import enum
import os
from os.path import dirname, abspath
import pathlib
from llvmlite import binding
import ast

load('helpers.sage')
load('boxspline.sage')
load('lattice.sage')
load('splinespace.sage')
load('codegen.sage')
load('horner.sage')


# Setup the LLVM context
import llvmlite.binding as llvm
llvm.initialize()
llvm.initialize_native_target()
llvm.initialize_native_asmprinter()

class EnumAction(argparse.Action):
    def __init__(self, **kwargs):
        enum_type = kwargs.pop("type", None)
        if enum_type is None:
            raise ValueError("type must be assigned an Enum when using EnumAction")
        if not issubclass(enum_type, enum.Enum):
            raise TypeError("type must be an Enum when using EnumAction")
        kwargs.setdefault("choices", tuple(e.value for e in enum_type))
        super(EnumAction, self).__init__(**kwargs)
        self._enum = enum_type

    def __call__(self, parser, namespace, values, option_string=None):
        value = self._enum(values)
        setattr(namespace, self.dest, value)


class FsActions(enum.Enum):
    DeriveSpline = "derive"
    DeriveVoronoi = "derive_voronoi"
    LatticeGen = "define_lattice"
    ComputeSplineSpace = "compute_spline_space"
    CodeGen = "codegen"
    RhoGen = "rho"


parser = argparse.ArgumentParser(description='Main fast-spline automation/experimentation CLI application')
parser.add_argument('action', metavar='action', type=FsActions, help=
    """The main verb/action for the script: Possible verbs are 
    [derive]: Derives a spline with --spline_matrix as the integer box spline matrix ex: "[-1,1,1;1,-1,1;1,1,-1] \n"
    [derive_voronoi]: Derives a voronoi spline with --spline one of  "vbcc2","vbcc3","vbcc4", "vfcc2", "vfcc3, "vfcc4" \n"
    """)
parser.add_argument('-sm', '--spline_matrix', metavar='spline_matrix', type=str, default=None,
                    help='This option is the input for the [derive_spline] verb, this should be a matrix specified in the form '
                         ' "[1,0,0,1; 0, 1, 1, 0; 0,0,1,1]". ')
parser.add_argument('-lm', '--lattice_matrix', metavar='lattice_matrix', type=str, default=None, help='lattice spline matrix')
parser.add_argument('-l', '--lattice', metavar='lattice', type=str, default=None, help='path to the lattice')
parser.add_argument('-s', '--spline', metavar='spline', type=str, default=None, help='path to the spline')
parser.add_argument('-ss', '--spline-space', metavar='spline', type=str, default=None, help='path to the spline space')
parser.add_argument('-o', '--output', metavar='output', type=str, default=None, help='path for output storage')
parser.add_argument('-r', '--rho', metavar='rho', type=str, default=None, help='path to the rho')

parser.add_argument('-rt', '--rho-type', metavar='rho-type', type=str, default=None, help='indicator or pp')
parser.add_argument('-rs', '--rho-shift', metavar='rho-shift', type=str, default=None,
                    help='s-dimensional vector defining the shift for rho')
parser.add_argument('-rpp', '--rho-pp', metavar='rho-pp', type=str, default=None,
                    help='matrix describing the parallelpiped for rho (if rho-type == pp)')

parser.add_argument('--quiet', metavar='', type=bool, default=False, help='omit console logging info')
parser.add_argument('--lite', metavar='',  type=bool, default=False, help='omit outputting human readable info to output directory')
parser.add_argument('--lg-disable-lattice-check', type=str, default=False, help='')
parser.add_argument('--ss-disable-reflective',  type=ast.literal_eval, default=False, help='')
parser.add_argument('--ss-shift-half', type=str, default='yes', help='')
parser.add_argument('--cg-disable-branchless',  type=str, default=False, help='')
parser.add_argument('--cg-strict-branchless',  type=ast.literal_eval, default=True, help='')
parser.add_argument('--cg-memfetch-planning', type=str, default=False, help='')
parser.add_argument('--cg-lanes', type=int, default=1, help='')
parser.add_argument('--cg-group-size',  type=int, default=1, help='')
parser.add_argument('--cg-pipeline-depth', metavar='cg_pipeline_depth',  type=int, default=1, help='')
parser.add_argument('--cg-fetches',  type=str, default="ABSTRACT", help='')
parser.add_argument('--cg-func-name',  type=str, default=None, help='')
parser.add_argument('--cg-target', type=str, default=binding.Target.from_default_triple().triple)
parser.add_argument('--cg-disable-reflective',  type=ast.literal_eval, default=False, help='')
parser.add_argument('--region-cache-path', type=str, default=None)

args = parser.parse_args()


def parse_matrix(matrix_str, require_int=False):
    num = lambda _: ZZ(int(_))
    rows = matrix_str.replace("[", "").replace("]", "").split(';')
    rows = [[num(p.strip()) for p in _.split(',')] for _ in rows]
    if len(rows) == 0:
        raise ValueError("Invalid Matrix!")
    if any([len(_) != len(rows[0]) for _ in rows]):
        raise ValueError("Invalid Matrix!")

    if require_int:
        return matrix(ZZ, rows)
    return matrix(rows)


def parse_shift(shift_str, require_int=False):
    num = lambda _: QQ(_)
    els = shift_str.replace("[", "").replace("]", "").split(',')
    els = [num(p.strip()) for p in els]
    # print(els)

    if len(els) == 0:
        raise ValueError("Invalid Vector!")

    if require_int:
        return vector(ZZ, els)
    return vector(els)

def create_dir(dir):
    if len(dir.split(os.path.pathsep)[-1].split('.')) > 1 and not dir[-1]==os.pathsep:
        dir =  dirname(abspath(dir))
    dir = abspath(dir)
    pathlib.Path(dir).mkdir(parents=True, exist_ok=True)
    return dir


def load_obj(dir):
    if os.path.isdir(dir):
        return load(os.path.join(dir, 'object'))
    return load(os.path.join(dirname(dir), 'object'))


def has_args(args, *args_):
    for _ in args_:
        if getattr(args, _) is None:
            print(f"Argument '{_}' required!")
            print(f"Invalid usage, use -h to display the help")
            return False
    return True


def read_vs(spline_in):
    x_0, x_1, x_2 = [var('x_%d' % _)  for _ in range(3)]
    vs = []
    with open(spline_in, 'r') as fin:
        for verts in fin:
            verts = sage_eval(verts)
            pp = sage_eval(fin.readline(), locals={'x_0': x_0, 'x_1': x_1, 'x_2': x_2})
            vs += [(Polyhedron(vertices=verts), pp)]
    return vs


def main(args):
    if args.action == FsActions.DeriveSpline:

        if not has_args(args, 'spline_matrix', 'output'):
            exit()
        try:
            matrix = parse_matrix(args.spline_matrix, True)
        except ValueError:
            raise ValueError(f"Spline matrix: '{args.spline_matrix}' is incorrectly formatted")

        output_dir = create_dir(args.output)
        print(f"Deriving spline:\n{matrix}\nOutput: {output_dir}\nInit...")
        bs = BoxSpline(matrix.columns(), verbose=True)

        print(f"Deriving PP regions...")
        pp_regions = bs.get_pp_regions()

        print(f"Ok, saving...")
        save(bs, os.path.join(output_dir, 'object'))

        with open(os.path.join(output_dir, 'object.json.spline'), 'w') as fout:
            fout.write(':)')

        # TODO: Print out human readible descriptions
        # TODO: Print out json parsable descriptions
    elif args.action == FsActions.RhoGen:
        matrix = None
        if args.rho_type is None:
            raise ValueError('You must specify the type of region for rho')

        if args.rho_type == "indicator":
            r_type = "indicator"
        elif args.rho_type == "pp":
            r_type = "pp"
            if args.rho_pp is None:
                raise ValueError("You must specify --rho-pp if you use --rho-type pp")

            try:
                matrix = parse_matrix(args.rho_pp, True)
            except ValueError:
                raise ValueError(f"Parallelpiped matrix: '{args.rho_pp}' is incorrectly formatted")
        else:
            raise ValueError(f"Invalid rho type '{args.rho_type}'")

        shift = None
        if args.rho_shift is not None:
            shift = parse_shift(args.rho_shift)
            if r_type == "pp":
                if len(shift) != len(matrix.columns()):
                    raise ValueError("dimension mistmatch between PP and shift vector!")

        output_dir = create_dir(args.output)
        save({"type": r_type, "matrix": matrix, "shift": shift}, os.path.join(output_dir, 'object'))
        with open(os.path.join(output_dir, 'object.json.rho'), 'w') as fout:
            fout.write(f"type: {r_type}\n")
            fout.write(f"matrix: {matrix}\n")
            fout.write(f"shift: {shift}\n")

    elif args.action == FsActions.DeriveVoronoi:
        print(f"voronoi spline.... for {args.spline}")
        output_dir = create_dir(args.output)

        if args.spline in ['bcc2', 'bcc3', 'bcc4', 'fcc2', 'fcc3', 'fcc4']:
            spline = read_vs('voronoi/%s.txt' % args.spline)
        else:
            raise ValueError()

        print("Opened spline, caching...")
        save(spline, os.path.join(output_dir, 'object'))
        with open(os.path.join(output_dir, 'object.json.spline'), 'w') as fout:
            fout.write(':)')

    elif args.action == FsActions.LatticeGen:
        if not has_args(args, 'lattice_matrix', 'output'):
            exit()
        try:
            matrix = parse_matrix(args.lattice_matrix, True)
        except ValueError:
            raise ValueError(f"lattice matrix: '{args.lattice_matrix}' is incorrectly formatted")

        output_dir = create_dir(args.output)
        print(f"Deriving lattice for gmat:\n{matrix}\nOutput: {output_dir}\nInit...")
        lattice = IntegerLattice(matrix)

        print(f"Ok, saving... {matrix.det()}")
        save(lattice, os.path.join(output_dir, 'object'))

        with open(os.path.join(output_dir, 'object.json.lattice'), 'w') as fout:
            fout.write(':)')

        # TODO: Print out human readible descriptions
        # TODO: Print out json parsable descriptions

    elif args.action == FsActions.ComputeSplineSpace:
        if not has_args(args, 'spline', 'lattice', 'output'):
            exit()
        spline = load_obj(args.spline)
        lattice = load_obj(args.lattice)
        rho = load_obj(args.rho)

        rho_shift = rho["shift"]
        rho = rho["type"] if rho["type"] == "indicator" else rho["matrix"]

        output_dir = create_dir(args.output)
        ss = SplineSpace(spline,
                         lattice,
                         enable_reflective=not args.ss_disable_reflective,
                         rho=rho,
                         rho_shift=rho_shift,
                         )
        save(ss, os.path.join(output_dir, 'object'))
        # TODO: Print out human readible descriptions
        # TODO: Print out json parsable descriptions
        with open(os.path.join(output_dir, 'object.json.ss'), 'w') as fout:
            fout.write(':)')

        if args.region_cache_path is not None:
            ecg = EvalCodeGen(ss,
                              lookup_type="ABSTRACT",
                              lane_count=args.cg_lanes,
                              group_size=args.cg_group_size,
                              pipeline_depth=args.cg_pipeline_depth,
                              target_triple=args.cg_target,

                              cache_path=args.region_cache_path)

            ecg = EvalCodeGen(ss,
                              lookup_type="CPU_LIN",
                              lane_count=args.cg_lanes,
                              group_size=args.cg_group_size,
                              pipeline_depth=args.cg_pipeline_depth,
                              target_triple=args.cg_target,
                              cache_path=args.region_cache_path)


    elif args.action == FsActions.CodeGen:

        if not has_args(args, 'spline_space', 'output'):
            exit()
        ss = load_obj(args.spline_space)
        output_dir = create_dir(args.output)

        ecg = EvalCodeGen(ss,
                          lookup_type=args.cg_fetches,
                          lane_count=args.cg_lanes,
                          group_size=args.cg_group_size,
                          pipeline_depth=args.cg_pipeline_depth,
                          target_triple=args.cg_target,
                          cache_path=args.region_cache_path,
                          strict_branchless=args.cg_strict_branchless,
                          reconstruction_name=args.cg_func_name)
        l = ecg.llvm_generate_vec()
        print('saving to ', os.path.join(output_dir, f"generated_{args.cg_target}.ll"))
        with open(os.path.join(output_dir, f"generated_{args.cg_target}.ll"), 'w') as output:
            output.write(l)


        # TODO: Print out json parsable summary


if __name__ == "__main__":
    main(args)
