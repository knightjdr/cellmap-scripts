const mockFS = require('mock-fs');
const fs = require('fs');

const writeArray = require('./write-array');

const mockedFileSystem = {};
mockFS(mockedFileSystem);

const expected = `term\tname\tic\tcount
a\tA\t0\t10
b\tB\t0.5\t5
c\tC\t1\t1
`;

afterAll(() => {
  mockFS.restore();
});

describe('Write an array to file', () => {
  it('should write an array of string to a file', async (done) => {
    const arr = [
      { name: 'A', term: 'a', ic: 0, count: 10 },
      { name: 'B', term: 'b', ic: 0.5, count: 5 },
      { name: 'C', term: 'c', ic:1, count: 1 },
    ];
    await writeArray(arr, './file.txt');
    const data = fs.readFileSync('./file.txt', 'utf8');
    expect(data).toBe(expected);
    done();
  });
});
