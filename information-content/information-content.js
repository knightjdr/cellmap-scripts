const round = require('./round');

const informationContent = (numGenes, frequency, termMap) => (
  Object.entries(frequency).map(([term, count]) => ({
    count,
    ic: round(-Math.log10(count / numGenes), 3),
    name: termMap[term] || 'unknown',
    term,
  }))
);

module.exports = informationContent;
