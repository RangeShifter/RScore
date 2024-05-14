#include "GeneticFitnessTrait.h"

// ----------------------------------------------------------------------------------------
//  Initialisation constructor
//  Called when initialising community
//  Sets up initial values, and immutable attributes (distributions and parameters)
//	that are defined at the species-level
// ----------------------------------------------------------------------------------------
GeneticFitnessTrait::GeneticFitnessTrait(SpeciesTrait* P)
{
	pSpeciesTrait = P;
	ExpressionType expressionType = pSpeciesTrait->getExpressionType();

	initialise();

	_inherit_func_ptr = (pSpeciesTrait->getPloidy() == 1) ? &GeneticFitnessTrait::inheritHaploid : &GeneticFitnessTrait::inheritDiploid; //this could be changed if we wanted some alternative form of inheritance

	DistributionType mutationDistribution = pSpeciesTrait->getMutationDistribution();
	map<GenParamType, float> mutationParameters = pSpeciesTrait->getMutationParameters();

	switch (mutationDistribution) {
	case UNIFORM:
	{
		if (mutationParameters.count(MAX) != 1)
			throw logic_error("Error:: genetic load mutation uniform distribution parameter must contain one max value (e.g. max= ) \n");
		if (mutationParameters.count(MIN) != 1)
			throw logic_error("Error:: genetic load mutation uniform distribution parameter must contain one min value (e.g. min= ) \n");
		break;
	}
	case NORMAL:
	{
		if (mutationParameters.count(MEAN) != 1)
			throw logic_error("Error:: genetic load mutation distribution set to normal so parameters must contain one mean value (e.g. mean= ) \n");
		if (mutationParameters.count(SD) != 1)
			throw logic_error("Error:: genetic load mutation distribution set to normal so parameters must contain one sdev value (e.g. sdev= ) \n");
		break;
	}
	case GAMMA:
	{
		if (mutationParameters.count(SHAPE) != 1)
			throw logic_error("Error:: genetic load mutation distribution set to gamma so parameters must contain one shape value (e.g. shape= ) \n");
		if (mutationParameters.count(SCALE) != 1)
			throw logic_error("Error:: genetic load mutation distribution set to gamma so parameters must contain one scale value (e.g. scale= ) \n");
		break;
	}
	case NEGEXP:
	{
		if (mutationParameters.count(MEAN) != 1)
			throw logic_error("Error:: genetic load mutation distribution set to negative exponential (negative decay) so parameters must contain one mean value (e.g. mean= ) \n");
		break;
	}
	default:
		throw logic_error("Error:: wrong parameter value for genetic load mutation model, must be uniform/normal/gamma/negExp \n");
	}

	DistributionType dominanceDistribution = pSpeciesTrait->getDominanceDistribution();
	map<GenParamType, float> dominanceParameters = pSpeciesTrait->getDominanceParameters();

	switch (dominanceDistribution) {
	case UNIFORM:
	{
		if (dominanceParameters.count(MAX) != 1)
			throw logic_error("Error:: genetic load dominance uniform distribution parameter must contain one max value (e.g. max= ) \n");
		if (dominanceParameters.count(MIN) != 1)
			throw logic_error("Error:: genetic load dominance uniform distribution parameter must contain one min value (e.g. min= ) \n");
		break;
	}
	case NORMAL:
	{
		if (dominanceParameters.count(MEAN) != 1)
			throw logic_error("Error:: genetic load dominance distribution set to normal so parameters must contain one mean value (e.g. mean= ) \n");
		if (dominanceParameters.count(SD) != 1)
			throw logic_error("Error:: genetic load dominance distribution set to normal so parameters must contain one sdev value (e.g. sdev= ) \n");
		break;
	}
	case GAMMA:
	{
		if (dominanceParameters.count(SHAPE) != 1)
			throw logic_error("Error:: genetic load dominance distribution set to gamma so parameters must contain one shape value (e.g. shape= ) \n");
		if (dominanceParameters.count(SCALE) != 1)
			throw logic_error("Error:: genetic load dominance distribution set to gamma so parameters must contain one scale value (e.g. scale= ) \n");
		break;
	}
	case NEGEXP:
	{
		if (dominanceParameters.count(MEAN) != 1)
			throw logic_error("Error:: genetic load dominance distribution set to negative exponential (negative decay) so parameters must contain mean value (e.g. mean= ) \n");
		break;
	}
	case SCALED:
	{
		break; // no parameter
	}
	default:
	{
		throw logic_error("Error:: wrong parameter value for genetic load dominance model, must be uniform/normal/gamma/negExp/scaled \n");
		break;
	}
	}
}

// ----------------------------------------------------------------------------------------
// Inheritance constructor
// Copies immutable features from a parent trait
// Only called via clone()
// ----------------------------------------------------------------------------------------
GeneticFitnessTrait::GeneticFitnessTrait(const GeneticFitnessTrait& T) : pSpeciesTrait(T.pSpeciesTrait), _inherit_func_ptr(T._inherit_func_ptr)
{}

void GeneticFitnessTrait::initialise() {
	// All positions start at wild type, mutations accumulate through simulation
	const set<int> genePositions = pSpeciesTrait->getGenePositions();
	short ploidy = pSpeciesTrait->getPloidy();
	const vector<shared_ptr<Allele>> wildTypeGene(ploidy, wildType);
	for (auto position : genePositions) {
		genes.insert(make_pair(position, wildTypeGene));
	}
}

// ----------------------------------------------------------------------------------------
// Mutate uniform
// ----------------------------------------------------------------------------------------
void GeneticFitnessTrait::mutate()
{
	const int positionsSize = pSpeciesTrait->getPositionsSize();
	const auto& genePositions = pSpeciesTrait->getGenePositions();
	const short ploidy = pSpeciesTrait->getPloidy();
	const float mutationRate = pSpeciesTrait->getMutationRate();
	float newSelectionCoef;
	float newDominanceCoef;

	auto rng = pRandom->getRNG();

	for (int p = 0; p < ploidy; p++) {
		// Determine nb of mutations
		unsigned int NbMut = pRandom->Binomial(positionsSize, mutationRate);

		if (NbMut > 0) {
			vector<int> mutationPositions;
			// Draw which positions mutate
			sample(genePositions.begin(), genePositions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {
				auto it = genes.find(m);
				if (it == genes.end())
					throw runtime_error("Locus sampled for mutation doesn't exist.");
				do {
					newSelectionCoef = drawSelectionCoef();
				} while (!pSpeciesTrait->isValidTraitVal(newSelectionCoef));
				do {
					newDominanceCoef = drawDominance(newSelectionCoef);
				} while (newDominanceCoef < 0.0);
				it->second[p] = make_shared<Allele>(newSelectionCoef, newDominanceCoef);
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// get dominance value for new mutation
// ----------------------------------------------------------------------------------------
float GeneticFitnessTrait::drawDominance(float selCoef) {

	DistributionType dominanceDistribution = pSpeciesTrait->getDominanceDistribution();
	map<GenParamType, float> dominanceParameters = pSpeciesTrait->getDominanceParameters();

	float h;
	switch (dominanceDistribution) {
	case UNIFORM:
	{
		float maxD = dominanceParameters.find(MAX)->second;
		float minD = dominanceParameters.find(MIN)->second;
		h = pRandom->FRandom(minD, maxD);
		break;
	}
	case NORMAL:
	{
		const float mean = dominanceParameters.find(MEAN)->second;
		const float sd = dominanceParameters.find(SD)->second;
		do {
			h = static_cast<float>(pRandom->Normal(mean, sd));
		} while (h <= 0.0);
		break;
	}
	case GAMMA:
	{
		const float shape = dominanceParameters.find(SHAPE)->second;
		const float scale = dominanceParameters.find(SCALE)->second;
		h = static_cast<float>(pRandom->Gamma(shape, scale));
		break;
	}
	case NEGEXP:
	{
		const float mean = dominanceParameters.find(MEAN)->second;
		h = static_cast<float>(pRandom->NegExp(mean));
		break;
	}
	case SCALED:
	{
		const float min = 0;
		const float max = static_cast<float>(exp((-log(2 * 0.36) / 0.05) * selCoef));
		h = static_cast<float>(pRandom->FRandom(min, max));
		break;
	}

	default:
	{
		throw logic_error("Error:: wrong parameter value for genetic load dominance model, must be uniform/normal/gamma/negExp/scaled \n");
		break;
	}
	}
	return h;
}

// ----------------------------------------------------------------------------------------
// Get selection coefficient for new mutation
// 
// Selection coefficients will usually be between 0 and 1, but may,
// if the mutation distribution enable it, take a negative value
// down to -1 representing the effect of beneficial mutations
// ----------------------------------------------------------------------------------------
float GeneticFitnessTrait::drawSelectionCoef() {

	DistributionType mutationDistribution = pSpeciesTrait->getMutationDistribution();
	map<GenParamType, float> mutationParameters = pSpeciesTrait->getMutationParameters();

	float s = 0.0; // default selection coefficient is 0

	switch (mutationDistribution) {
	case UNIFORM:
	{
		float maxD = mutationParameters.find(MAX)->second;
		float minD = mutationParameters.find(MIN)->second;
		s = pRandom->FRandom(minD, maxD); // no check here, min and max should already be constrained to valid values
		break;
	}
	case NORMAL:
	{
		const float mean = mutationParameters.find(MEAN)->second;
		const float sd = mutationParameters.find(SD)->second;
		do {
			s = static_cast<float>(pRandom->Normal(mean, sd));
		} while (!pSpeciesTrait->isValidTraitVal(s));
		break;
	}
	case GAMMA:
	{
		const float shape = mutationParameters.find(SHAPE)->second;
		const float scale = mutationParameters.find(SCALE)->second;
		do {
			s = static_cast<float>(pRandom->Gamma(shape, scale));
		} while (!pSpeciesTrait->isValidTraitVal(s));
		break;
	}
	case NEGEXP:
	{
		const float mean = mutationParameters.find(MEAN)->second;
		do {
			s = static_cast<float>(pRandom->NegExp(mean));
		} while (!pSpeciesTrait->isValidTraitVal(s));
		break;
	}
	default:
	{
		throw logic_error("Error:: wrong parameter value for genetic load mutation model, must be uniform/normal/gamma/negExp/scaled \n");
		break;
	}
	}
	return s;
}


// ----------------------------------------------------------------------------------------
//  Wrapper to inheritance function
// ----------------------------------------------------------------------------------------
void GeneticFitnessTrait::inheritGenes(const bool& fromMother, TTrait* parentTrait, set<unsigned int> const& recomPositions, int startingChromosome)
{
	auto parentCast = dynamic_cast<GeneticFitnessTrait*> (parentTrait); // must convert TTrait to GeneticFitnessTrait
	const auto& parent_seq = parentCast->getGenes();
	(this->*_inherit_func_ptr) (fromMother, parent_seq, recomPositions, startingChromosome);
}

// ----------------------------------------------------------------------------------------
// Inheritance for diploid, sexual species
// Called once for each parent. Given a list of recombinant sites, 
// populates offspring genes with appropriate parent alleles
// Assumes mother genes are inherited first
// ----------------------------------------------------------------------------------------
void GeneticFitnessTrait::inheritDiploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome) {

	auto it = recomPositions.lower_bound(parentGenes.begin()->first);
	int nextBreakpoint = *it;
	auto distance = std::distance(recomPositions.begin(), it);
	if (distance % 2 != 0)
		parentChromosome = 1 - parentChromosome; // switch chromosome

	for (auto const& [locus, allelePair] : parentGenes) {

		// Switch chromosome if locus is past recombination site
		while (locus > nextBreakpoint) {
			std::advance(it, 1);
			nextBreakpoint = *it;
			parentChromosome = 1 - parentChromosome;
		}

		if (locus <= nextBreakpoint) {
			auto& parentAllele = allelePair[parentChromosome];
			auto itGene = genes.find(locus);
			if (itGene == genes.end()) {
				// locus does not exist yet, create and initialise it
				if (!fromMother) throw runtime_error("Father-inherited locus does not exist.");
				vector<shared_ptr<Allele>> newAllelePair(2); // always diploid
				newAllelePair[sex_t::FEM] = parentAllele;
				genes.insert(make_pair(locus, newAllelePair));
			} 
			else { // father, locus already exists
				if (fromMother) throw runtime_error("Mother-inherited locus already exists.");
				itGene->second[sex_t::MAL] = parentAllele;
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// Inheritance for haploid, asexual species
// Simply pass down parent genes
// Arguments are still needed to match overloaded function in base class
// ----------------------------------------------------------------------------------------
void GeneticFitnessTrait::inheritHaploid(const bool& fromMother, map<int, vector<shared_ptr<Allele>>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome)
{
	genes = parentGenes;
}

// ----------------------------------------------------------------------------------------
// Expression genetic load
// ----------------------------------------------------------------------------------------
float GeneticFitnessTrait::express() {

	float phenotype = 1.0;  // base chance of viability
	float sA, sB, hA, hB, sumDomCoeffs, hLocus;

	for (auto const& [locus, pAllelePair] : genes)
	{
		shared_ptr<Allele> pAlleleA  = pAllelePair[0] == 0 ? wildType : pAllelePair[0];
		shared_ptr<Allele> pAlleleB = pAllelePair[1] == 0 ? wildType : pAllelePair[1];

		sA = pAlleleA->getAlleleValue();
		hA = pAlleleA->getDominanceCoef();
		if (pSpeciesTrait->getPloidy() == 2) {
			sB = pAlleleB->getAlleleValue();
			hB = pAlleleB->getDominanceCoef();
		}
		else {
			sB = 0.0;
			hB = 0.0;
		}
		sumDomCoeffs = hA + hB;
		hLocus = sumDomCoeffs == 0.0 ? 0.0 : hA / sumDomCoeffs;
		phenotype *= 1 - hLocus * sA - (1 - hLocus) * sB;
	}
	return phenotype;
}

// ----------------------------------------------------------------------------------------
// Check if specific locus is heterozygote
// ----------------------------------------------------------------------------------------
bool GeneticFitnessTrait::isHeterozygoteAtLocus(int locus) const {
	// assumes diploidy
	auto it = genes.find(locus);
	if (it == genes.end())
		throw runtime_error("Genetic load gene queried for heterozygosity does not exist.");
	else {
		shared_ptr<Allele> alleleRight = it->second[0] == 0 ? wildType : it->second[0];
		shared_ptr<Allele> alleleLeft = it->second[1] == 0 ? wildType : it->second[1];
		return alleleRight != alleleLeft;
	}
}

// ----------------------------------------------------------------------------------------
// Count heterozygote loci in genome 
// ----------------------------------------------------------------------------------------
int GeneticFitnessTrait::countHeterozygoteLoci() const {
	// assumes diploidy
	int count = 0;
	for (auto const& [locus, allelePair] : genes) {
		shared_ptr<Allele> alleleLeft = allelePair[0] == 0 ? wildType : allelePair[0];
		shared_ptr<Allele> alleleRight = allelePair[1] == 0 ? wildType : allelePair[1];
		count += alleleLeft != alleleRight;
	}
	return count;
}

// ----------------------------------------------------------------------------------------
// Get allele value at locus
// ----------------------------------------------------------------------------------------
float GeneticFitnessTrait::getAlleleValueAtLocus(short whichChromosome, int position) const {

	auto it = genes.find(position);
	if (it == genes.end())
		throw runtime_error("The genetic load locus queried for its allele value does not exist.");
	return it->second[whichChromosome] == 0 ? wildType->getAlleleValue() : it->second[whichChromosome]->getAlleleValue();
}