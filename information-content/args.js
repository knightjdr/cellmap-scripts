const args = () => {
  const options = {
    annotations: '',
    namespace: 'C',
    obo: '',
  };
  process.argv.forEach((arg) => {
    if (arg.includes('--annotations')) {
      options.annotations = arg.split('=')[1];
    } else if (arg.includes('--namespace')) {
      options.namespace = arg.split('=')[1];
    } else if (arg.includes('--obo')) {
      options.obo = arg.split('=')[1];
    }
  });
  return options;
};

module.exports = args;
