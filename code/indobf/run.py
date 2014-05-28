#!/usr/bin/env sage -python

from __future__ import print_function

from bp_sww import SWWBranchingProgram
from bp_barrington import BarringtonBranchingProgram
from test import test_circuit

import argparse, os, sys, time

def bp(args):
    testdir = 'circuits'
    if args.barrington:
        bpclass = BarringtonBranchingProgram
    else:
        bpclass = SWWBranchingProgram
    if args.test_circuit:
        test_circuit(args.test_circuit, bpclass, None, False, args)
    if args.test_all:
        for circuit in os.listdir('circuits'):
            path = os.path.join(testdir, circuit)
            if os.path.isfile(path) and path.endswith('.circ'):
                test_circuit(path, bpclass, None, False, args)
    if args.load_circuit:
        bp = bpclass(args.load_circuit, verbose=args.verbose)
        if args.obliviate:
            bp.obliviate()
        if args.eval:
            r = bp.evaluate(args.eval)
            print('Output = %d' % r)

def obf(args):
    from obfuscator import BarringtonObfuscator, SWWObfuscator
    if args.barrington:
        bpclass = BarringtonBranchingProgram
        obfclass = BarringtonObfuscator
    else:
        bpclass = SWWBranchingProgram
        obfclass = SWWObfuscator
    if args.test_circuit:
        test_circuit(args.test_circuit, bpclass, obfclass, True, args)
    else:
        directory = None
        if args.load_obf:
            directory = args.load_obf
        elif args.load_circuit:
            start = time.time()
            obf = obfclass(verbose=args.verbose)
            directory = args.save if args.save \
                        else '%s.obf.%d' % (args.load_circuit, args.secparam)
            obf.obfuscate(args.load_circuit, args.secparam, directory,
                          obliviate=args.obliviate)
            end = time.time()
            print("Obfuscation took: %f seconds" % (end - start))
            obf.cleanup()
        else:
            print("One of --load-obf, --load-circuit, or --test-circuit must be used")
            sys.exit(1)

        if args.eval:
            assert directory
            obf = obfclass(verbose=args.verbose)
            r = obf.evaluate(directory, args.eval)
            print('Output = %d' % r)

def main():
    parser = argparse.ArgumentParser(
        description='Run indistinguishability obfuscator.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()

    parser_bp = subparsers.add_parser(
        'bp',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help='commands for circuit -> branching program conversion')
    parser_bp.add_argument('--eval', metavar='INPUT', type=str, action='store',
                           help='evaluate branching program on INPUT')
    parser_bp.add_argument('--load-circuit', metavar='FILE', type=str,
                           action='store', help='load circuit from FILE')
    parser_bp.add_argument('--test-circuit', metavar='FILE', type=str,
                           action='store',
                           help='test FILE circuit -> bp conversion')
    parser_bp.add_argument('--test-all', action='store_true',
                           help='test circuit -> bp conversion')
    parser_bp.add_argument('--secparam', metavar='N', type=int, action='store',
                           default=24, help='security parameter')

    parser_bp.add_argument('--obliviate', action='store_true',
                           help='obliviate the branching program')
    parser_bp.add_argument('--barrington', action='store_true',
                           help="use Barrington's approach")
    parser_bp.add_argument('-v', '--verbose', action='store_true',
                           help='be verbose')
    parser_bp.set_defaults(func=bp)

    parser_obf = subparsers.add_parser(
        'obf',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help='commands for obfuscating a circuit/branching program')
    parser_obf.add_argument('--eval', metavar='INPUT', type=str, action='store',
                            help='evaluate obfuscation on INPUT')
    parser_obf.add_argument('--load-obf', metavar='DIR', type=str,
                            action='store',
                            help='load obfuscation from DIR')
    parser_obf.add_argument('--load-circuit', metavar='FILE', type=str,
                            action='store',
                            help='load circuit from FILE and obfuscate')
    parser_obf.add_argument('--test-circuit', metavar='FILE', type=str,
                            action='store',
                            help='test circuit from FILE')
    parser_obf.add_argument('--save', metavar='DIR', type=str, action='store',
                            help='save obfuscation to DIR')
    parser_obf.add_argument('--secparam', metavar='N', type=int, action='store',
                            default=24, help='security parameter')

    parser_obf.add_argument('--obliviate', action='store_true',
                            help='obliviate the branching program')
    parser_obf.add_argument('--barrington', action='store_true',
                            help="use Barrington's approach")
    parser_obf.add_argument('-v', '--verbose', action='store_true', 
                            help='be verbose')
    parser_obf.set_defaults(func=obf)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass
