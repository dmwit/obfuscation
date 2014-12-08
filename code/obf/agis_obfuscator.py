import _obfuscator as _obf
from obfuscator import Obfuscator
from agis_bp import AGISBranchingProgram

import os, time

def to_long(lst):
    return [long(i) for i in lst]

def pad(array, length, bplength):
    if len(array) < length:
        zeros = [to_long([0] * bplength)]
        return array + (zeros * (length - len(array)))
    else:
        return array

class AGISObfuscator(Obfuscator):
    def __init__(self, mlm='CLT', verbose=False):
        super(AGISObfuscator, self).__init__(_obf, mlm=mlm, verbose=verbose)
        self._mlm_dict = {
            'CLT': self._gen_clt_mlm_params,
            'GGH': self._gen_ggh_mlm_params,
        }

    def _gen_clt_mlm_params(self, secparam, kappa, width, nzs, directory):
        self.logger('Generating MLM parameters...')
        start = time.time()
        if not os.path.exists(directory):
            os.mkdir(directory)
        self._state, primes = _obf.setup_clt(secparam, kappa, width, nzs,
                                             directory)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return primes

    def _gen_ggh_mlm_params(self, secparam, kappa, width, nzs, directory):
        self.logger('Generating GGH parameters...')
        start = time.time()
        self._state, prime = _obf.setup_ggh(secparam, kappa, width, directory)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return [prime]

    def _gen_mlm_params(self, secparam, kappa, width, nzs, directory):
        return self._mlm_dict[self._mlm](secparam, kappa, width, nzs, directory)

    def _construct_bps(self, bpclass, num, circuit, obliviate):
        self.logger('Constructing %d BP%s...' % (num, '' if num == 1 else 's'))
        start = time.time()
        bps = []
        for _ in xrange(num):
            bp = bpclass(circuit, verbose=self._verbose, obliviate=obliviate)
            bp.set_straddling_sets()
            bps.append(bp)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return bps

    def _construct_bookend_vectors(self, bps, num, nzs):
        def compute_vectors(bps, num):
            start = time.time()
            bplength = len(bps[0].s)
            ss = pad([to_long(bp.s * bp.m0i) for bp in bps], num, bplength)
            ts = pad([to_long(bp.m0 * bp.t) for bp in bps], num, bplength)
            end = time.time()
            self.logger('  Computing bookend vectors: %f' % (end - start))
            return ss, ts
        self.logger('Constructing bookend vectors...')
        start = time.time()
        ss, ts = compute_vectors(bps, num)
        _obf.encode_vectors(self._state, ss, [nzs - 2], 's_enc')
        _obf.encode_vectors(self._state, ts, [nzs - 1], 't_enc')
        end = time.time()
        self.logger('Took: %f' % (end - start))

    def _randomize(self, secparam, bps, primes):
        self.logger('Randomizing BPs...')
        start = time.time()
        for bp, prime in zip(bps, primes):
            bp.randomize(prime, mlm=self._mlm)
            bp.set_straddling_sets()
        end = time.time()
        self.logger('Took: %f' % (end - start))

    def _obfuscate(self, bps, num):
        for i in xrange(len(bps[0])):
            self.logger('Obfuscating layer...')
            start = time.time()
            bplength = len(bps[0][i].zero.transpose().list())
            zeros = pad([to_long(bp[i].zero.transpose().list()) for bp in bps],
                        num, bplength)
            ones = pad([to_long(bp[i].one.transpose().list()) for bp in bps],
                       num, bplength)
            nrows = bps[0][i].zero.nrows()
            ncols = bps[0][i].zero.ncols()
            _obf.encode_layers(self._state, i, nrows, ncols, bps[0][i].inp,
                               zeros, ones, bps[0][i].zeroset, bps[0][i].oneset)
            end = time.time()
            self.logger('Took: %f' % (end - start))

    def obfuscate(self, circuit, secparam, directory, obliviate=False,
                  nslots=None, kappa=None):
        start = time.time()

        self._remove_old(directory)
        if nslots is None:
            nslots = secparam

        # create a dummy branching program to determine parameters
        bp = AGISBranchingProgram(circuit, verbose=self._verbose,
                                  obliviate=obliviate)
        # add two to kappa due to the bookend vectors
        if not kappa:
            kappa = len(bp) + 2
        # construct straddling sets, and add two to the number of Zs to take
        # bookend vectors into account
        nzs = bp.set_straddling_sets() + 2
        # width is the column/row-length of the matrices
        width = bp.size

        primes = self._gen_mlm_params(secparam, kappa, width, nzs, directory)
        num = min(nslots, len(primes))

        bps = self._construct_bps(AGISBranchingProgram, num, circuit, obliviate)
        self._randomize(secparam, bps, primes)
        self._construct_bookend_vectors(bps, num, nzs)
        self._obfuscate(bps, num)

        end = time.time()
        self.logger('Obfuscation took: %f' % (end - start))
        if self._verbose:
            _obf.max_mem_usage()

    def evaluate(self, directory, inp):
        return self._evaluate(directory, inp, _obf.evaluate_agis, _obf)

    def cleanup(self):
        pass
