#ifdef UNIT_TESTS

// Unit tests for the weighted dispersal kernel (perception accuracy, rejection sampling)

#include "../Individual.h"
#include "../Landscape.h"
#include "../Species.h"
#include "../Cell.h"
#include <cmath>
#include <limits>

// Species with one stage, no stage/sex dependence, constant transfer mortality 0
static void setUpKernelSpeciesWK(Species& sp, const float meanDist, const int resol) {
	demogrParams d{};
	d.stageStruct = false;
	sp.setDemogr(d);
	transferRules trfr{}; // value-initialised: all flags false
	trfr.usesMovtProc = false;
	trfr.indVar = trfr.sexDep = trfr.stgDep = false;
	trfr.twinKern = trfr.distMort = false;
	sp.setTrfrRules(trfr);
	sp.setFullKernel(false);
	trfrMortParams mort{};
	mort.fixedMort = 0.0;
	sp.setMortParams(mort);
	settleRules sett{};
	sett.wait = false;
	sp.setSettRules(0, 0, sett);
	trfrKernelParams kern;
	kern.meanDist1 = meanDist;
	kern.meanDist2 = meanDist;
	kern.probKern1 = 1.0;
	sp.setSpKernTraits(0, 0, kern, resol);
}

static bool throwsInvalidArgument(void (*f)(void)) {
	try { f(); }
	catch (invalid_argument&) { return true; }
	return false;
}

// A1: parameters, defaults and range checks
void testPerceptionParameters() {
	Species sp;
	assert(sp.getSpPerception(0, 0).perceptionAcc == 0.0f);
	assert(!sp.usesPerception());
	trfrPerceptionRules r = sp.getPerceptionRules();
	assert(r.maxRejections == 100);
	assert(!r.transferFile);
	assert(std::isinf(sp.getPerceptionDist()));
	assert(sp.getTransferRef() == 1.0f);

	trfrPerceptionParams p;
	p.perceptionAcc = 0.7f;
	sp.setSpPerception(0, 1, p);
	assert(sp.getSpPerception(0, 1).perceptionAcc == 0.7f);
	assert(sp.getSpPerception(0, 0).perceptionAcc == 0.0f);
	assert(sp.usesPerception());
	p.perceptionAcc = 1.0f; // upper bound is allowed
	sp.setSpPerception(0, 0, p);
	assert(sp.getSpPerception(0, 0).perceptionAcc == 1.0f);
	// invalid stage/sex: ignored, getter returns the default
	sp.setSpPerception(-1, 0, p);
	assert(sp.getSpPerception(gMaxNbStages, 0).perceptionAcc == 0.0f);

	assert(throwsInvalidArgument([] { Species s; trfrPerceptionParams q; q.perceptionAcc = -0.1f; s.setSpPerception(0, 0, q); }));
	assert(throwsInvalidArgument([] { Species s; trfrPerceptionParams q; q.perceptionAcc = 1.1f; s.setSpPerception(0, 0, q); }));
	assert(throwsInvalidArgument([] { Species s; trfrPerceptionParams q; q.perceptionAcc = std::numeric_limits<float>::quiet_NaN(); s.setSpPerception(0, 0, q); }));
	assert(throwsInvalidArgument([] { Species s; trfrPerceptionRules q; q.maxRejections = 0; s.setPerceptionRules(q); }));
	assert(throwsInvalidArgument([] { Species s; trfrPerceptionRules q; q.transferFile = true; s.setPerceptionRules(q); }));
	r.maxRejections = 20;
	sp.setPerceptionRules(r);
	assert(sp.getPerceptionRules().maxRejections == 20);

	sp.overridePerceptionExtension(250.0f, 0.4f);
	assert(sp.getPerceptionDist() == 250.0f);
	assert(sp.getTransferRef() == 0.4f);
}

// A2: transfer values derived from the landscape
void testTransferValues() {
	Species sp;
	sp.createHabK(3);
	sp.setHabK(0, 0.0);
	sp.setHabK(1, 50.0);
	sp.setHabK(2, 100.0);
	assert(sp.getHabTransfer(0) == 0.0f);
	assert(sp.getHabTransfer(1) == 0.5f);
	assert(sp.getHabTransfer(2) == 1.0f);
	assert(sp.getHabTransfer(-1) == 0.0f);
	assert(sp.getHabTransfer(3) == 0.0f);

	Cell cHab(0, 0, nullptr, 1);
	assert(sp.getTransferValue(&cHab, 0, 0, 3) == 0.5f);
	Cell cQual(0, 0, nullptr, 40.0f); // habitat quality 40 %
	assert(std::abs(sp.getTransferValue(&cQual, 0, 2, 1) - 0.4f) < 1e-6f);
	// % cover: one layer, 60 % cover of a habitat with K = max K
	Species sp1;
	sp1.createHabK(1);
	sp1.setHabK(0, 80.0);
	Cell cCover(0, 0, nullptr, 60.0f);
	assert(std::abs(sp1.getTransferValue(&cCover, 0, 1, 1) - 0.6f) < 1e-6f);

	assert(sp.getTransferMax() == 1.0f); // default before any landscape: safe bound
}

// A3: T_max of the landscape for w_max
void testTransferMax() {
	// habitat codes: T is normalised by max K, so T_max = 1 (0 if no habitat has K > 0)
	{
		landParams ls_params = createDefaultLandParams(3);
		Landscape ls;
		ls.setLandParams(ls_params, true);
		ls.setCellArray();
		ls.addCellToLand(new Cell(0, 0, nullptr, 0));
		ls.addCellToLand(new Cell(1, 1, nullptr, 1));
		Species sp;
		sp.createHabK(2);
		sp.setHabK(0, 0.0);
		sp.setHabK(1, 20.0);
		sp.updateTransferMax(&ls, 0);
		assert(sp.getTransferMax() == 1.0f);
		sp.setHabK(1, 0.0);
		sp.updateTransferMax(&ls, 0);
		assert(sp.getTransferMax() == 0.0f);
	}
	// habitat quality: highest quality in the landscape / 100
	{
		landParams ls_params = createDefaultLandParams(3);
		ls_params.rasterType = 2;
		ls_params.nHab = ls_params.nHabMax = 1;
		Landscape ls;
		ls.setLandParams(ls_params, true);
		ls.setCellArray();
		ls.addCellToLand(new Cell(0, 0, nullptr, 20.0f));
		ls.addCellToLand(new Cell(2, 1, nullptr, 60.0f));
		ls.addCellToLand(new Cell(1, 2, nullptr, 35.0f));
		assert(ls.getMaxHabitatQuality(0) == 60.0f);
		Species sp;
		sp.createHabK(1);
		sp.setHabK(0, 100.0);
		sp.updateTransferMax(&ls, 0);
		assert(std::abs(sp.getTransferMax() - 0.6f) < 1e-6f);
	}
}

// A3: effective attempt cap N_eff = min(ceil(pi * (ln 4 * meandist)^2), MaxRejections)
void testAttemptCap() {
	Species sp;
	assert(sp.kernelAttemptCap(1.0) == 7);   // delta = 100 m, 100 m cells
	assert(sp.kernelAttemptCap(4.0) == 97);  // delta = 100 m, 25 m cells
	assert(sp.kernelAttemptCap(20.0) == 100); // limited by MaxRejections
	assert(sp.kernelAttemptCap(0.01) == 1);  // at least 1
	trfrPerceptionRules r;
	r.maxRejections = 50;
	sp.setPerceptionRules(r);
	assert(sp.kernelAttemptCap(4.0) == 50);
}

// A4: direction of dispersal in a landscape with an unsuitable west half (T = 0)
// and a suitable east half (T = 1); returns the share of landings in the east half
static double eastShare(const float aS, long& nWest, long& nFallback, const int nInds) {
	const int dim = 21, centre = 10;
	landParams ls_params = createDefaultLandParams(dim);
	Landscape ls;
	ls.setLandParams(ls_params, true);
	ls.setCellArray();
	for (int x = 0; x < dim; x++)
		for (int y = 0; y < dim; y++)
			ls.addCellToLand(new Cell(x, y, nullptr, x >= centre ? 1 : 0));

	Species sp;
	sp.createHabK(2);
	sp.setHabK(0, 0.0);
	sp.setHabK(1, 100.0);
	setUpKernelSpeciesWK(sp, 3.0f, ls_params.resol);
	trfrPerceptionParams p;
	p.perceptionAcc = aS;
	sp.setSpPerception(0, 0, p);
	ls.allocatePatches(&sp);
	ls.updateCarryingCapacity(&sp, 0, 0);
	sp.updateTransferMax(&ls, 0);

	Cell* start = ls.findCell(centre, centre);
	Patch* startPatch = start->getPatch();
	sp.resetKernelFallbacks();
	long nEast = 0;
	nWest = 0;
	for (int i = 0; i < nInds; i++) {
		Individual ind(&sp, start, startPatch, 1, 0, 0, 0.0, false, 0);
		ind.moveKernel(&ls, &sp, false, 0);
		Cell* c = ind.getCurrCell();
		assert(c != nullptr); // reflecting boundaries: every individual lands
		if (c->getLocn().x >= centre) nEast++; else nWest++;
	}
	nFallback = sp.getKernelFallbacks();
	return (double)nEast / nInds;
}

void testWeightedKernelDirection() {
	const int n = 3000;
	long nWest, nFallback;
	const double share0 = eastShare(0.0f, nWest, nFallback, n);
	assert(nFallback == 0); // unweighted: no attempts, no fallback
	assert(share0 > 0.4 && share0 < 0.7); // roughly symmetric (centre column counts as east)

	const double shareHalf = eastShare(0.5f, nWest, nFallback, n);
	const double share1 = eastShare(1.0f, nWest, nFallback, n);
	// with a_s = 1, unsuitable cells (T = 0) are never accepted: west landings only by fallback
	assert(nWest == nFallback);
	assert(share1 > 0.95);
	// a_s = 0.5 lies in between
	assert(shareHalf > share0 + 0.05 && shareHalf < share1);
}

// A4: fallback when no proposal can be accepted, and at the 1000-step safety cap
void testWeightedKernelFallback() {
	// (1) only the natal cell is suitable; a_s = 1: every proposal has w = 0
	{
		const int dim = 21, centre = 10;
		landParams ls_params = createDefaultLandParams(dim);
		Landscape ls;
		ls.setLandParams(ls_params, true);
		ls.setCellArray();
		for (int x = 0; x < dim; x++)
			for (int y = 0; y < dim; y++)
				ls.addCellToLand(new Cell(x, y, nullptr, (x == centre && y == centre) ? 1 : 0));
		Species sp;
		sp.createHabK(2);
		sp.setHabK(0, 0.0);
		sp.setHabK(1, 100.0);
		setUpKernelSpeciesWK(sp, 3.0f, ls_params.resol);
		trfrPerceptionParams p;
		p.perceptionAcc = 1.0f;
		sp.setSpPerception(0, 0, p);
		ls.allocatePatches(&sp);
		ls.updateCarryingCapacity(&sp, 0, 0);
		sp.updateTransferMax(&ls, 0);
		Cell* start = ls.findCell(centre, centre);
		sp.resetKernelFallbacks();
		const int n = 200;
		for (int i = 0; i < n; i++) {
			Individual ind(&sp, start, start->getPatch(), 1, 0, 0, 0.0, false, 0);
			ind.moveKernel(&ls, &sp, false, 0);
			assert(ind.getCurrCell() != nullptr && ind.getCurrCell() != start);
			assert(ind.getStatus() == 6); // lands in the matrix and dies, as without weighting
		}
		assert(sp.getKernelFallbacks() == n); // every jump ends in the fallback
	}
	// (2) safety cap: all cells suitable but with T = 0.01 (K = 1 against a max K of 100 in one
	// corner cell), a_s = 1, large MaxRejections and a kernel mean much larger than the landscape,
	// so that most proposals fall outside (reflecting boundaries) and the 1000-step cap is
	// reached; the last valid proposal must then be used instead of status 6
	{
		const int dim = 11, centre = 5;
		landParams ls_params = createDefaultLandParams(dim);
		Landscape ls;
		ls.setLandParams(ls_params, true);
		ls.setCellArray();
		for (int x = 0; x < dim; x++)
			for (int y = 0; y < dim; y++)
				ls.addCellToLand(new Cell(x, y, nullptr, (x == 0 && y == 0) ? 2 : 1));
		Species sp;
		sp.createHabK(3);
		sp.setHabK(0, 0.0);
		sp.setHabK(1, 1.0);
		sp.setHabK(2, 100.0);
		setUpKernelSpeciesWK(sp, 200.0f, ls_params.resol);
		trfrPerceptionParams p;
		p.perceptionAcc = 1.0f;
		sp.setSpPerception(0, 0, p);
		trfrPerceptionRules r;
		r.maxRejections = 100000;
		sp.setPerceptionRules(r);
		ls.allocatePatches(&sp);
		ls.updateCarryingCapacity(&sp, 0, 0);
		sp.updateTransferMax(&ls, 0);
		Cell* start = ls.findCell(centre, centre);
		sp.resetKernelFallbacks();
		const int n = 50;
		for (int i = 0; i < n; i++) {
			Individual ind(&sp, start, start->getPatch(), 1, 0, 0, 0.0, false, 0);
			ind.moveKernel(&ls, &sp, false, 0);
			assert(ind.getCurrCell() != nullptr);
			assert(ind.getStatus() == 2); // potential settler, not status 6
		}
		assert(sp.getKernelFallbacks() > 0);
	}
}

void testWeightedKernel() {
	testPerceptionParameters();
	testTransferValues();
	testTransferMax();
	testAttemptCap();
	testWeightedKernelDirection();
	testWeightedKernelFallback();
}

#endif // UNIT_TESTS
