const tiers = (numGenes) => {
  console.log('Tiers:');
  console.log(`tier 1 - 2% ${Math.ceil(numGenes * 0.02)}`);
  console.log(`tier 2 - 10% ${Math.ceil(numGenes * 0.1)}`);
  console.log(`tier 3 - 25% ${Math.ceil(numGenes * 0.25)}`);
  console.log(`tier 4 - 100% ${Math.ceil(numGenes)}`);
};

module.exports = tiers;
