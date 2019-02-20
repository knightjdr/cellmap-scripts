const informationContent = require('./information-content');

describe('Information content', () => {
  it('should return an array of terms with information content', () => {
    const freq = {
      a: 2,
      b: 3,
      c: 1,
    };
    const map = {
      a: 'A',
      b: 'B',
      c: 'C',
    };
    const expected = [
      { ic: 0.301, name: 'A', term: 'a' },
      { ic: 0.125, name: 'B', term: 'b' },
      { ic: 0.602, name: 'C', term: 'c' },
    ];
    expect(informationContent(4, freq, map)).toEqual(expected);
  });
});
