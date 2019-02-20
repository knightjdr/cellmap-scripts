const fs = require('fs');

const writeArray = (arr, file) => (
  new Promise((resolve) => {
    const stream = fs.createWriteStream(file);

    stream.write('term\tname\tic\tcount\n');
    arr.forEach((item) => {
      stream.write(`${item.term}\t${item.name}\t${item.ic}\t${item.count}\n`);
    });
    stream.end();

    stream.on('finish', () => {
      resolve();
    });
  })
);

module.exports = writeArray;
