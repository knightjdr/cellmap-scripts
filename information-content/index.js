const args = require('./args');
const informationContent = require('./information-content');
const tiers = require('./tiers');
const writeArray = require('./write-array');
const { readAnnotations } = require('./read-annotations');
const { readObo } = require('./read-obo');

const options = args();

let map;

readObo(options.obo)
  .then((obo) => {
    map = obo.map;
    return readAnnotations(options.annotations, options.namespace, obo.parents);
  })
  .then(([numGenes, frequency]) => {
    const ic = informationContent(numGenes, frequency, map);
    tiers(numGenes);
    writeArray(ic, './information-content.txt');
  })
  .catch((err) => {
    console.log(err);
  });
