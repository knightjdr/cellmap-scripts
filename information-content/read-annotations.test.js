const mockFS = require('mock-fs');

const { addParents, parseGenes, termFrequency } = require('./read-annotations');

const goAnnotations = `UniProtKB	Q9BUL8	PDCD10		GO:0000139	GO_REF:0000039	IEA	UniPr	C
UniProtKB	Q9BUL8	PDCD10		GO:0001525	GO_REF:0000037	IEA	UniProtKB-KW:KW-0037	P
UniProtKB	Q9BUL8	PDCD10		GO:0005515	PMID:16189514	IPI	UniProtKB:O00506	C
UniProtKB	Q9BUL8	PDCD10		GO:0005515	PMID:16189514	IPI	UniProtKB:Q9Y6E0	C
UniProtKB	Q9P289	STK26		GO:0000139	GO_REF:0000107	IEA	UniProtKB:F1LXV3|ensembl:ENSR	C
UniProtKB	Q9P289	STK26		GO:0000287	PMID:11641781	IDA		F
UniProtKB	Q9P289	STK26		GO:0004672	PMID:11641781	IDA		F
UniProtKB	Q9P289	STK26		GO:0004672	PMID:22797597	IMP		F
UniProtKB	Q9P289	STK26		GO:0005515	PMID:15037601	IPI	UniProtKB:Q08379	C
UniProtKB	Q9P289	STK26		GO:0005515	PMID:17360971	IPI	UniProtKB:Q9BUL8	C
`;

const mockedFileSystem = {
  'annotations.txt': goAnnotations,
};
mockFS(mockedFileSystem);

afterAll(() => {
  mockFS.restore();
});

describe('Add parents', () => {
  it('should all parent terms to list of GO terms', () => {
    const annotations = {
      PDCD10: ['GO:0000139', 'GO:0001525', 'GO:0005515'],
      STK26: ['GO:0000139', 'GO:0000287', 'GO:0004672', 'GO:0005515'],
    };
    const parents = {
      'GO:0000139': ['GO:AAAAAAA', 'GO:BBBBBBB'],
      'GO:0001525': ['GO:AAAAAAA', 'GO:CCCCCCC', 'GO:DDDDDDD'],
      'GO:0005515': ['GO:BBBBBBB', 'GO:CCCCCCC'],
    };
    const expected = {
      PDCD10: [
        'GO:0000139',
        'GO:AAAAAAA',
        'GO:BBBBBBB',
        'GO:0001525',
        'GO:CCCCCCC',
        'GO:DDDDDDD',
        'GO:0005515',
      ],
      STK26: [
        'GO:0000139',
        'GO:AAAAAAA',
        'GO:BBBBBBB',
        'GO:0000287',
        'GO:0004672',
        'GO:0005515',
        'GO:CCCCCCC',
      ],
    };
    expect(addParents(annotations, parents)).toEqual(expected);
  });
});

describe('Parse genes', () => {
  it('should get all GO terms from a namespace in an annotation file, but not duplicates', () => {
    const expected = {
      PDCD10: ['GO:0000139', 'GO:0005515'],
      STK26: ['GO:0000139', 'GO:0005515'],
    };
    return expect(parseGenes('annotations.txt', 'C')).resolves.toEqual(expected);
  });
});

describe('Term frequency', () => {
  let results;

  beforeAll(() => {
    const annotations = {
      1: ['a'],
      2: ['a', 'b', 'c'],
      3: ['c' , 'd'],
    };
    results = termFrequency(annotations);
  });

  it('should sum the number of occurrences of a value in an object of arrays', () => {
    const expected = {
      a: 2,
      b: 1,
      c: 2,
      d: 1,
    };
    expect(results[1]).toEqual(expected);
  });

  it('should return the number of genes in the dataset', () => {
    expect(results[0]).toEqual(3);
  });
});
