/* 序列工坊 app shell service worker. Requires HTTPS (or localhost). */
const CACHE_PREFIX = 'seq-workshop-next-' + encodeURIComponent(self.registration.scope) + '-';
const CACHE_NAME = CACHE_PREFIX + 'v4-57935b17164c6695';
const INDEX_URL = new URL('./index.html', self.registration.scope).href;
const SHELL_URLS = [
  new URL('./', self.registration.scope).href,
  INDEX_URL,
  new URL('./manifest.json', self.registration.scope).href,
  new URL('./assets/echarts.min.js', self.registration.scope).href,
  new URL('./assets/vendor/three.min.js', self.registration.scope).href,
  new URL('./assets/structures/structures.js', self.registration.scope).href,
  new URL('./assets/js/structure3d.js', self.registration.scope).href,
  new URL('./assets/js/sequence-core.js', self.registration.scope).href,
  new URL('./assets/js/workspace-format.js', self.registration.scope).href,
  new URL('./assets/js/product-workspace.js', self.registration.scope).href,
  new URL('./assets/js/app-bundle.js', self.registration.scope).href,
  new URL('./assets/js/module-config.js', self.registration.scope).href,
  new URL('./assets/js/module-router.js', self.registration.scope).href,
  new URL('./assets/css/base.css', self.registration.scope).href,
  new URL('./assets/css/module-shell.css', self.registration.scope).href,
  new URL('./icons/icon.svg', self.registration.scope).href,
  new URL('./icons/icon-192.png', self.registration.scope).href,
  new URL('./icons/icon-512.png', self.registration.scope).href
];

self.addEventListener('install', event => {
  event.waitUntil(
    caches.open(CACHE_NAME)
      .then(cache => cache.addAll(SHELL_URLS))
      .then(() => self.skipWaiting())
  );
});

self.addEventListener('activate', event => {
  event.waitUntil(
    caches.keys()
      .then(keys => Promise.all(
        keys
          .filter(key => key.startsWith(CACHE_PREFIX) && key !== CACHE_NAME)
          .map(key => caches.delete(key))
      ))
      .then(() => self.clients.claim())
  );
});

self.addEventListener('message', event => {
  if (event.data && event.data.type === 'SKIP_WAITING') {
    self.skipWaiting();
  }
});

self.addEventListener('fetch', event => {
  const request = event.request;
  if (request.method !== 'GET') return;

  const url = new URL(request.url);
  // The app has no API or runtime CDN dependency; only same-origin app assets are cached.
  if (url.origin !== self.location.origin) return;

  if (request.mode === 'navigate') {
    event.respondWith(networkFirstNavigation(request));
    return;
  }

  event.respondWith(cacheFirstAsset(request));
});

async function networkFirstNavigation(request) {
  try {
    const response = await fetch(request);
    if (response.ok) {
      const cache = await caches.open(CACHE_NAME);
      await cache.put(INDEX_URL, response.clone());
    }
    return response;
  } catch (error) {
    const cache = await caches.open(CACHE_NAME);
    const cached = await cache.match(request) || await cache.match(INDEX_URL);
    if (cached) return cached;
    return new Response(
      '序列工坊暂时无法加载，请恢复网络后重试。',
      {status: 503, headers: {'Content-Type': 'text/plain; charset=utf-8'}}
    );
  }
}

async function cacheFirstAsset(request) {
  const cache = await caches.open(CACHE_NAME);
  const cached = await cache.match(request);
  if (cached) return cached;

  const response = await fetch(request);
  if (response.ok && response.type === 'basic') {
    await cache.put(request, response.clone());
  }
  return response;
}
