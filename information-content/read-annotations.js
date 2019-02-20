const LineByLineReader = require('line-by-line');

const arrayUnique = require('./array-unique');

// Read genes and their annotations from GO annotations file (go-basic.ogo)
const parseGenes = (file, namespace) => (
  new Promise((resolve, reject) => {
    const genes = {};

    const lineReader = new LineByLineReader(file);
    lineReader.on('line', (line) => {
      const [,, gene, , term, , , , ns] = line.split('\t');
      if (ns === namespace) {
        if (genes[gene]) {
          if (!genes[gene].includes(term)) {
            genes[gene].push(term);
          }
        } else {
          genes[gene] = [term];
        }
      }
    });
    lineReader.on('end', () => {
      resolve(genes);
    });
    lineReader.on('error', (err) => {
      reject(err);
    });
  })
);

// Add all parent terms to a genes's annotation list
const addParents = (annotations, parents) => {
  const fullAnnotations = annotations;
  Object.entries(annotations).forEach(([gene, terms]) => {
    fullAnnotations[gene] = [];
    terms.forEach((term) => {
      const parentTerms = parents[term] || [];
      fullAnnotations[gene] = [
        ...fullAnnotations[gene],
        term,
        ...parentTerms,
      ];
    });
    fullAnnotations[gene] = arrayUnique(fullAnnotations[gene]);
  });
  return fullAnnotations;
};

// Calculate how often a GO term occurs based on an object of genes with annotations
// and also return the number of genes in the dataset.
const termFrequency = (annotations) => {
  const termFrequency = {};
  Object.values(annotations).forEach((terms) => {
    terms.forEach((term) => {
      termFrequency[term] = termFrequency[term] ? termFrequency[term] + 1 : 1;
    });
  });
  return [Object.keys(annotations).length, termFrequency];
};

const readAnnotations = (file, namespace, parents) => (
  new Promise((resolve, reject) => {
    parseGenes(file, namespace)
      .then((annotations) => {
        const annotationWithParents = addParents(annotations, parents);
        resolve(termFrequency(annotationWithParents));
      })
      .catch((err) => {
        reject(err);
      });
  })
);

module.exports = {
  addParents,
  parseGenes,
  readAnnotations,
  termFrequency,
};
